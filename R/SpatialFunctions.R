# Functions for working with spatially-resolved sequencing data

######################################################################################################
#######                                                                                         ######
#######                                  User-Facing Functions                                  ######
#######                                                                                         ######
######################################################################################################


#' Build polygons from a CSV file of boundary points created with defineKeyAreas
#' 
#' @param    path_to_polygons   Path to CSV file containing polygon boundary points
#' @param    coords             Data frame of coordinates with columns loc_name, x, and y
#'
#' @return   data frame with columns x, y, and loc_name representing all points within the defined polygons
#' @export
#' 
buildPolygons <- function(path_to_polygons, coords){

  # --- Check Parameters --- #
  if (!file.exists(path_to_polygons)){
    stop("path_to_polygons must be a valid file path")
  }
  if (!is.data.frame(coords) | !(all(c("loc_name", "x", "y") %in% colnames(coords)))){
    stop("coords must be a data frame with columns loc_name, x, and y")
  }

  # --- Read in polygon boundary points --- #
  polygons_boundary <- read.csv(path_to_polygons)
  polygons_boundary$x <- round(polygons_boundary$x, 3)
  polygons_boundary$y <- round(polygons_boundary$y, 3)
  polygons_boundary_coords <- polygons_boundary %>% select(x, y)

  # --- Separate polygons --- #
  # We may have multiple polygons, we need to look for them by repetitions of coordinates (since coords shouldn't be repeated except to close a polygon)
  polygons_boundaries <- list()
  n_left <- nrow(polygons_boundary_coords)
  temp_coords <- polygons_boundary_coords
  while(n_left > 0) {
    first_coord_idx <- 1
    repeat_coord <- which(temp_coords$x == temp_coords$x[first_coord_idx] & temp_coords$y == temp_coords$y[first_coord_idx])[2]
    if (length(repeat_coord) == 0) {
      stop("No repeat coordinates found, please check the polygon boundary coordinates. The first and last points of a polygon must be identical.")
    }
    polygons_boundaries[[length(polygons_boundaries) + 1]] <- temp_coords[1:repeat_coord, ]
    temp_coords <- temp_coords[-(1:repeat_coord), ]
    n_left <- nrow(temp_coords)
  }
  
  # --- Encapsulate all internal points in polygons --- #
  polygons <- list()
  for (i in 1:length(polygons_boundaries)) {
    polygons[[i]] <- define_polygon(spatial_grid = coords,
                                      by_all_or_boundary = "boundary",
                                      by_names_or_coords = "coords",
                                      names_or_coords = polygons_boundaries[[i]])
    polygons[[i]]$is_boundary <- ifelse(polygons[[i]]$loc_name %in% polygons_boundary$ID[polygons_boundary$polygon == i], TRUE, FALSE)
  }

  # Collapse list and number polygons
  for (i in 1:length(polygons)) {
    polygons[[i]]$polygon_id <- i
  }
  polygons <- dplyr::bind_rows(polygons)
  return(polygons)
}


#' Create centroids for spatial downsampling using modified farthest point sampling algorithm to exclude boundaries
#' 
#' @param    spatial_grid                          Data frame with columns x, y, and loc_name covering the spatial grid
#' @param    exclude_points                        Data frame with columns x, y, and loc_name covering all points to exclude for building centroids
#' @param    num_centroids                         Number of centroids to downsample non-key area points onto (default 100)
#' @param    num_neighbors_non_boundary            Number of neighboring points needed to be considered non-boundary
#' @param    neighbor_distance_threshold           Distance threshold to consider points as neighbors
#' @param    seed                                  Seed for centroid placement with FPS (default 123456)
#' @param    plot_excluded_points_and_boundaries   Logical indicating whether to plot key areas and boundaries for visual inspection
#' @param    point_size                            Size of points in the plot if plot_excluded_points_and_boundaries is TRUE
#'
#' @return   data frame with columns x, y, and loc_name representing all points within the defined polygons
#' @export
#'
createCentroidsFPS <- function(spatial_grid, exclude_points = NULL, num_centroids = 100, num_neighbors_non_boundary, neighbor_distance_threshold, seed = 123456, plot_excluded_points_and_boundaries = FALSE, point_size = 2) {

  # --- Check Parameters --- #
  if (!is.data.frame(spatial_grid) | !(all(c("x", "y", "loc_name") %in% colnames(spatial_grid)))){
    stop("spatial_grid must be a data frame with columns x, y, and loc_name")
  }
  spatial_grid <- spatial_grid %>% select(x, y, loc_name)
  
  if (!is.null(exclude_points)){
    if (!is.data.frame(exclude_points) | !(all(c("x", "y", "loc_name") %in% colnames(exclude_points)))){
      stop("exclude_points must be a data frame with columns x, y, and loc_name")
    }
    exclude_points <- exclude_points %>% select(x, y, loc_name)
  }

  if (!is.numeric(num_centroids) | num_centroids <= 0 | num_centroids %% 1 != 0){
    stop("num_centroids must be a positive integer")
  }

  if (!is.numeric(num_neighbors_non_boundary) | num_neighbors_non_boundary <= 0 | num_neighbors_non_boundary %% 1 != 0){
    stop("num_neighbors_non_boundary must be a positive integer")
  }

  if (!is.numeric(neighbor_distance_threshold) | neighbor_distance_threshold <= 0){
    stop("neighbor_distance_threshold must be a positive real number")
  }

  # --- Build Centroids with Modified FPS --- #
  # Remove any points to exclude from the spatial grid
  coords <- spatial_grid
  if (!is.null(exclude_points)) {
    coords <- coords %>% dplyr::anti_join(exclude_points, by = c("x", "y"))
  }
  coords <- coords %>% dplyr::select(x, y)
  
  # Identify boundary points
  coords_dist <- stats::dist(coords) %>% as.matrix()
  coords$name <- rownames(coords) 
  coords$boundary <- ifelse(apply(coords_dist, 1, function(x) sum(x < neighbor_distance_threshold) < (num_neighbors_non_boundary + 1)), 1, 0) # + 1 to account for the point itself
  boundary_points <- coords %>% dplyr::filter(boundary == 1) %>% dplyr::select(x, y)
  
  # Plot the excluded points and boundaries if requested
  if (plot_excluded_points_and_boundaries){
    g <- ggplot2::ggplot() + ggplot2::geom_point(data = spatial_grid, aes(x=x, y=y), color = "black", size = 4) + ggplot2::theme_bw() +
           ggplot2::theme(panel.grid = element_blank()) +
           ggplot2::geom_point(data = boundary_points, aes(x = x, y = y), color = "purple", size = point_size) +
           ggplot2::theme(legend.position = "none")
    if (!is.null(exclude_points)){
      g <- g + ggplot2::geom_point(data = exclude_points, aes(x = x, y = y), color = "green2", size = point_size)
    }
    print(g)
  }

  # Filter out boundary points before FPS
  interior_points <- spatial_grid %>% dplyr::anti_join(boundary_points, by = c("x", "y"))
  if (!is.null(exclude_points)) {
    interior_points <- interior_points %>% dplyr::anti_join(exclude_points, by = c("x", "y"))
  }
  interior_points <- interior_points %>% dplyr::select(x, y) %>% as.matrix()
  
  # Make sure we're not downsampling to more locations than it's possible to downsample to
  if (num_centroids > nrow(interior_points)) {
    stop("N (currently ", num_centroids, ") must be an integer less than or equal to the number of interior points of spatial_grid excluding the exclude_points.\n",
           "Number of interior points: ", nrow(interior_points))
  }
  
  # Run farthest point sampling on interior points
  set.seed(seed)
  initial_point_idx <- sample(1:nrow(interior_points), 1)
  downsampled_indices <- rdist::farthest_point_sampling(interior_points, metric = "euclidean", k = num_centroids, initial_point_index = initial_point_idx)

  # Retrieve the downsampled points and the loc_names
  downsampled_coords <- interior_points[downsampled_indices, ] %>% as.data.frame()
  downsampled_coords <- dplyr::left_join(downsampled_coords, spatial_grid, by = c("x", "y"))
  downsampled_coords$centroid_name <- paste0("centroid_", seq_len(nrow(downsampled_coords)))

  # Return as a data frame
  return(downsampled_coords)
  
}


#' Build a map from all points to the centroid they downsample onto via presence in Voronoi polygons drawn around centroids
#' 
#' @param    spatial_grid                          Data frame with columns x, y, and loc_name covering the spatial grid
#' @param    exclude_points                        Data frame with columns x, y, and loc_name covering all points to exclude for building centroids
#' @param    centroids                             Data frame with columns x, y, loc_name, and centroid_name covering the downsampled centroids
#' @param    num_neighbors_non_boundary            Number of neighboring points needed to be considered non-boundary
#' @param    rw_padding                            Padding to add to the rectangular window for Voronoi tessellation
#' @param    plot_centroids_voronoi_polygons       Logical indicating whether to plot centroids and Voronoi polygons for visual inspection
#' @param    point_size                            Size of points in the plot if plot_centroids_voronoi_polygons is TRUE
#'
#' @return   data frame mapping all points onto their respective downsampled centroids with columns loc.x, loc.y, loc_name, downsampled_centroid, centroid.x, centroid.y, centroid_loc_name, and area
#' @export
#'
buildDownsampleMap <- function(spatial_grid, exclude_points = NULL, centroids, rw_padding, plot_centroids_voronoi_polygons = FALSE, point_size = 2){

  # --- Check Parameters --- #
  if (!is.data.frame(spatial_grid) | !(all(c("x", "y", "loc_name") %in% colnames(spatial_grid)))){
    stop("spatial_grid must be a data frame with columns x, y, and loc_name")
  }
  spatial_grid <- spatial_grid %>% select(x, y, loc_name)
  
  if (!is.null(exclude_points)){
    if (!is.data.frame(exclude_points) | !(all(c("x", "y", "loc_name") %in% colnames(exclude_points)))){
      stop("exclude_points must be a data frame with columns x, y, and loc_name")
    }
    exclude_points <- exclude_points %>% select(x, y, loc_name)
    # Ensure exclude_points are a subset of spatial_grid
    if (!all(do.call(paste0, exclude_points) %in% do.call(paste0, spatial_grid))){
      stop("exclude_points must be a subset of spatial_grid")
    }
  }
    
  if (!is.data.frame(centroids) | !(all(c("x", "y", "loc_name", "centroid_name") %in% colnames(centroids)))){
    stop("centroids must be a data frame with columns x, y, loc_name, and centroid_name")
  }
  centroids <- centroids %>% select(x, y, loc_name, centroid_name)
  # Ensure that centroids are a subset of spatial_grid
  if (!all(do.call(paste0, centroids[, c("x", "y", "loc_name")]) %in% do.call(paste0, spatial_grid))){
    stop("centroids must all be in spatial_grid")
  }
  # Ensure that centroids are not in exclude points (if not null)
  if (!is.null(exclude_points)){
    if (any(do.call(paste0, centroids[, c("x", "y", "loc_name")]) %in% do.call(paste0, exclude_points))){
      stop("centroids cannot be in exclude_points")
    }
  }

  if (!is.numeric(rw_padding) | rw_padding <= 0){
    stop("rw_padding must be a positive real number for the Voronoi tessellation to work properly")
  }

  # --- Build Downsample Map --- #
  # Convert downsampled centroids to sf object
  centroids_sf <- sf::st_as_sf(data.frame(x = centroids$x, 
                                          y = centroids$y, 
                                          loc_name = centroids$loc_name,
                                          centroid_name = centroids$centroid_name),
                                 coords = c("x", "y"))

  # Compute Voronoi tessellation and extract the polygons
  rectangular_window <- c(min(spatial_grid$x) - rw_padding,  # should be xmin, xmax, ymin, ymax
                          max(spatial_grid$x) + rw_padding,
                          min(spatial_grid$y) - rw_padding,
                          max(spatial_grid$y) + rw_padding)
  voronoi_data <- deldir::deldir(centroids$x, centroids$y, rw = rectangular_window)
  voronoi_tiles <- deldir::tile.list(voronoi_data)

  # Convert to sf polygons
  voronoi_sf <- sf::st_sfc(
    lapply(voronoi_tiles, function(tile) {
      coords <- cbind(tile$x, tile$y)
      coords <- rbind(coords, coords[1, ])  # Ensure the polygon is closed
      sf::st_polygon(list(coords))
    })
  )

  # Assign the polygons to the centroids
  voronoi_sf <- sf::st_sf(geometry = voronoi_sf, data = centroids)

  # Plot the polygons and centroids if requested
  if (plot_centroids_voronoi_polygons){
    g <- ggplot2::ggplot() +
           ggplot2::geom_sf(data = voronoi_sf, fill = "transparent", color = "black") +
           ggplot2::geom_point(data = spatial_grid, aes(x = x, y = y), color = "black", size = point_size) +
           ggplot2::geom_point(data = centroids, aes(x = x, y = y), color = "red", size = point_size) +
           ggplot2::theme_bw() + ggplot2::theme(panel.grid = element_blank())
    if (!is.null(exclude_points)){
      g <- g + ggplot2::geom_point(data = exclude_points, aes(x = x, y = y), color = "green2", size = point_size)
    }
    print(g)
  }

  # Map points to downsample to centroids #
  # This is done by finding the centroid that each point is in the polygon of
  # In the case of ties, the point is assigned to the polygon with the fewest points currently in it
  # If this is still tied, we assign the point to the polygon with the smallest area
  # This is done to best balance the number of points in each polygon
  points_to_downsample <- spatial_grid 
  if (!is.null(exclude_points)){
      points_to_downsample <- points_to_downsample %>% dplyr::anti_join(exclude_points, by = c("x", "y"))
  }

  # Convert to an sf object
  points_to_downsample_sf <- sf::st_as_sf(points_to_downsample, coords = c("x", "y"))

  # Find which polygon each point is in
  polygon_idx <- sf::st_within(points_to_downsample_sf, voronoi_sf, sparse = FALSE) 
  polygon_idx <- matrix(as.numeric(polygon_idx), nrow = nrow(polygon_idx), ncol = ncol(polygon_idx))
  rownames(polygon_idx) <- points_to_downsample$loc_name; colnames(polygon_idx) <- centroids$centroid_name
  
  # Fill in the mapping data frame
  mapping <- points_to_downsample
  mapping$loc_name <- as.character(mapping$loc_name)
  mapping$downsampled_centroid <- sapply(mapping$loc_name, function(name) find_one(polygon_idx, name))

  # Take care of the points that are on the boundaries of polygons - have NAs in downsampled_centroid and intersect the boundaries of the polygons
  unassigned_points <- mapping %>% dplyr::filter(is.na(downsampled_centroid))
  if (nrow(unassigned_points) > 0){
    unassigned_points_sf <- sf::st_as_sf(unassigned_points, coords = c("x", "y"))
    unassigned_points_intersections <- sf::st_intersects(unassigned_points_sf, voronoi_sf, sparse = FALSE)
    centroid_table <- table(mapping$downsampled_centroid)
    
    # Determine which polygons the unassigned points intersect and assign them based on the rules above
    for (i in 1:nrow(unassigned_points)){
      unassigned_points_intersections <- sf::st_intersects(unassigned_points_sf[i, ], voronoi_sf, sparse = FALSE)
      potential_centroids <- centroids$centroid_name[which(unassigned_points_intersections == 1)]
      centroid_table_tmp <- centroid_table[potential_centroids]
      lowest_count_centroids <- names(centroid_table_tmp[which(centroid_table_tmp == min(centroid_table_tmp))])
      if (length(lowest_count_centroids) == 1){
        mapping$downsampled_centroid[mapping$loc_name == unassigned_points$loc_name[i]] <- lowest_count_centroids
      } else {
        # If there are multiple centroids with the same number of points, assign to the one with the smallest area
        areas <- sf::st_area(voronoi_sf[which(centroids$centroid_name %in% lowest_count_centroids), ])
        mapping$downsampled_centroid[mapping$loc_name == unassigned_points$loc_name[i]] <- lowest_count_centroids[which.min(areas)]
      }
    }
  }
  
  # Add the centroid coordinates to the mapping dataframe
  mapping <- dplyr::left_join(mapping, centroids, by = c("downsampled_centroid" = "centroid_name"))
  colnames(mapping) <- c("loc.x", "loc.y", "loc_name", "downsampled_centroid_name", "centroid.x", "centroid.y", "centroid_loc_name")
  mapping$area <- ifelse(mapping$loc_name %in% centroids$loc_name, "centroid", "downsampled")
  
  if (!is.null(exclude_points)){
    mapping_exclude_points <- data.frame(
          loc.x = exclude_points$x,
          loc.y = exclude_points$y,
          loc_name = as.character(exclude_points$loc_name),
          downsampled_centroid_name = paste0("exclude_point", 1:nrow(exclude_points)),
          centroid.x = exclude_points$x,
          centroid.y = exclude_points$y,
          centroid_loc_name = exclude_points$loc_name,
          area = "excluded_from_downsampling"
      )
      mapping <- dplyr::bind_rows(mapping, mapping_exclude_points)
  }
  mapping$area <- factor(mapping$area, levels = c("centroid", "downsampled", "excluded_from_downsampling"))
  
  return(mapping)
}


#' Converts spatial coordinates between pixel and cartesian coordinate systems
#' 
#' @param    coords                   data frame of coordinates
#' @param    x                        column in coords corresponding to x coordinate
#' @param    y                        column in coords corresponding to y coordinate
#' @param    direction                either "p2c" (pixel to cartesian) or "c2p" (cartesian to pixel)
#'
#' @return   data frame of converted coordinates with columns x and y
#' 
#' @export
convertCoordinates <- function(coords, x="imagerow", y="imagecol", direction="p2c") {

  # --- Check Parameters --- #
  if (!is.data.frame(coords)){
    stop("coords must be a data frame")
  }

  if (!all(c(x, y) %in% colnames(coords))){
    stop(paste0(x, " and ", y, " must be columns in coords"))
  }

  if (!direction %in% c("p2c", "c2p")){
    stop("direction must be either 'p2c' or 'c2p' for pixel to cartesian conversion or cartesian to pixel conversion")
  }

  # --- Convert Coordinates --- #
  # Remove other columns
  new_coords <- coords %>% select(-c(x, y))

  if (direction == "p2c") {  # Convert from pixel to cartesian
    new_coords$x <- coords[[x]]
    new_coords$y <- max(coords[[y]]) - coords[[y]]
  } else if (direction == "c2p") {  # Convert from cartesian to pixel
    new_coords$x <- coords[[x]]
    new_coords$y <- max(coords[[y]]) - coords[[y]]
  }
  return(new_coords)
}


#' Manually define key areas on a spatial gene expression plot
#' 
#' @param    seurat_obj         Spatial Seurat object containing coordinates and gene expression data
#' @param    assay              Assay containing gene expression data (default Spatial)
#' @param    dataset            Dataset within the assay containing gene expression data (default scale.data)
#' @param    autoload           Automatically build polygons and return them after closing the Shiny app (default TRUE)
#' @param    quiet              Logical indicating whether to suppress messages during Shiny app execution (default TRUE)
#'
#' @return   If autoload = TRUE, data frame containing all points (x, y, loc_name) within the key area polygon(s). Otherwise NULL
#' 
#' @export
defineKeyAreas <- function(seurat_obj, assay = "Spatial", dataset = "scale.data", autoload = TRUE, quiet = TRUE) {
  
  # --- Check Parameters --- #
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }

  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste0("Assay '", assay, "' not found in seurat_obj@assays."))
  }
  
  if (!dataset %in% slotNames(seurat_obj[[assay]])) {
    stop(paste0("Dataset '", dataset, "' not found in assay '", assay, "'."))
  }

  coords <- Seurat::GetTissueCoordinates(seurat_obj)

  # --- Run Shiny App --- #
  ui <- set_up_ui()
  
  server <- function(input, output, session) {
    if (quiet){
        old_warn <- getOption("warn")
        options(warn = -1)
    }

    seurat_obj <- seurat_obj
    
    # Track points and their polygon group
    selected_points <- shiny::reactiveVal(data.frame(cell = character(), x = numeric(), y = numeric(), polygon = integer(), stringsAsFactors = FALSE))

    # Track the current polygon number
    current_polygon <- shiny::reactiveVal(1)

    plot_data <- shiny::reactive({
      shiny::req(input$gene)
      plot_gene_expression_cartesian(seurat_obj, input$gene, assay = assay, dataset = dataset, point_size = input$point_size, option = "turbo")
    })
    
    output$gene_plot <- plotly::renderPlotly({
      base_plot <- plot_data()$plot
      selected_df <- selected_points()
      
      if (nrow(selected_df) > 0) {
        suppressWarnings({
          base_plot <- base_plot %>%
            plotly::add_markers(
              data = selected_df,
              x = ~x, y = ~y,
              marker = list(size = input$point_size * 5, color = 'white', line = list(width = 3, color = 'black')),
              name = "Selected Points"
            )
          
          # Draw lines within each polygon separately
          polygon_ids <- unique(selected_df$polygon)
          for (polygon in polygon_ids) {
            polygon_df <- selected_df[selected_df$polygon == polygon, ]
            if (nrow(polygon_df) > 1) {
              base_plot <- base_plot %>%
                plotly::add_trace(
                  data = polygon_df,
                  x = ~x, y = ~y,
                  mode = "lines",
                  line = list(color = "black", width = 2),
                  name = paste("Key Area", polygon),
                  showlegend = FALSE
                )
            }
          }
        })
      }
      base_plot
    })
    
    # Capture clicks and add closest point with current polygon ID
    shiny::observeEvent(plotly::event_data("plotly_click", source = "plot_click"), {
      click_data <- plotly::event_data("plotly_click", source = "plot_click")
      if (is.null(click_data)) return()
      
      expression_data <- plot_data()$data
      closest_point <- expression_data %>%
        dplyr::mutate(distance = sqrt((x - click_data$x)^2 + (y - click_data$y)^2)) %>%
        dplyr::arrange(distance) %>%
        dplyr::slice(1) %>%
        dplyr::select(-distance)
  
      # Add current polygon ID to point
      closest_point$polygon <- current_polygon()
      
      selected_points(rbind(selected_points(), closest_point))
    })
    
    # Show selected points in table (including polygon ID)
    output$selected_points <- shiny::renderTable({
      selected_points() %>% suppressWarnings(dplyr::select(-expression)) %>% dplyr::mutate(polygon = as.integer(polygon))
    })
    
    shiny::observeEvent(input$reset_points, {
      selected_points(data.frame(cell = character(), x = numeric(), y = numeric(), polygon = integer(), stringsAsFactors = FALSE))
      current_polygon(1)  # Reset polygon number too
    })
  
    shiny::observeEvent(input$remove_last_point, {
      current <- selected_points()
      if (nrow(current) > 0) {
        selected_points(current[-nrow(current), ])
      }
    })
    
    # Start a new polygon (increment the polygon counter)
    shiny::observeEvent(input$new_polygon, {
      current_polygon(current_polygon() + 1)
    })
  
    shiny::observeEvent(input$save_points, {
      write.csv(selected_points(), input$save_name, row.names = FALSE)
      shiny::showNotification("Polygons saved successfully!", type = "message")
      # Save and exit the app
      shiny::stopApp(returnValue = input$save_name)
      if (quiet){
          options(warn = old_warn)
      }
    })
  }
  
  options(shiny.launch.browser = TRUE)

  # Create the app object
  app <- shiny::shinyApp(ui, server)
  
  # Run the app and capture the file name
  filename <- shiny::runApp(app, launch.browser = TRUE)

  # --- After app closes, read in saved points and build polygons --- #
  print(paste("Saved polygons to", filename))

  if (autoload){
    print("Now building polygons...")

    coords_cartesian <- convertCoordinates(coords, x="imagecol", y="imagerow", direction="p2c")
    coords_cartesian$loc_name <- rownames(coords_cartesian)
    polygons <- buildPolygons(path_to_polygons = filename,
                                coords = coords_cartesian)
    return(polygons)
  } else {
    return(NULL)
  }
}


#' Manually define areas to remove on a spatial gene expression plot. Wrapper for defineKeyAreas for clarity
#' 
#' @param    seurat_obj         Spatial Seurat object containing coordinates and gene expression data
#' @param    assay              Assay containing gene expression data (default Spatial)
#' @param    dataset            Dataset within the assay containing gene expression data (default scale.data)
#'
#' @return   data frame containing all points (x, y, loc_name) within the remove area polygon(s)
#' 
#' @export
#' 
# Just an alias to defineKeyAreas for removing areas
defineRemoveAreas <- function(seurat_obj, assay = "Spatial", dataset = "scale.data") {
  defineKeyAreas(seurat_obj, assay = assay, dataset = dataset)
}

#' Downsample a clone assignment according to a downsample map
#' 
#' @param    clones             A tibble of \code{airrClone} objects, the output of \link{formatClones}
#' @param    downsample_map     Data frame mapping all points onto their respective downsampled centroids, the output of \link{buildDownsampleMap}
#' @param    method             Either "coords" to match by x and y coordinates, or "locs" to match by location name (default "coords")
#' @param    x                  If method is "coords", the name of the column in \code{f_clone} corresponding to x coordinates
#' @param    y                  If method is "coords", the name of the column in \code{f_clone} corresponding to y coordinates
#' @param    loc_name           If method is "locs", the name of the column in \code{f_clone} corresponding to location names
#' 
#' @return   A tibble of \code{airrClone} objects with locations downsampled according to the downsample map
#' 
#' @export
#' 
downsampleSpaceClones <- function(clones, 
                              downsample_map,
                              method = "coords",
                              x = NULL,
                              y = NULL,
                              loc_name = NULL,
                              ...) {

  # --- Check Parameters --- #
  if (!tibble::is_tibble(f_clones) | !all(sapply(f_clones$data, function(clone) inherits(clone, "airrClone")))){
    stop("f_clone must be a tibble of airrClone objects")
  }
  
  if (!is.data.frame(downsample_map) | !(all(c("loc.x", "loc.y", "loc_name", "downsampled_centroid_name", "centroid.x", "centroid.y", "centroid_loc_name", "area", "key_area_name") %in% colnames(downsample_map)))){
    stop("downsample_map must be a data frame with columns loc.x, loc.y, loc_name, downsampled_centroid_name, centroid.x, centroid.y, centroid_loc_name, area, and key_area_name" )
  }

  if (!method %in% c("coords", "locs")){
    stop("method must be either 'coords' or 'locs'")
  }

  if (method == "coords"){
    if (is.null(x) | is.null(y)){
      stop("If method is 'coords', x and y must be specified")
    }
    if (!all(sapply(clones$data, function(clone) all(c(x, y) %in% colnames(clone@data))))){
      stop(paste0(x, " and ", y, " must be columns in clones' data"))
    }
  } else if (method == "locs"){
    if (is.null(loc_name)){
      stop("If method is 'locs', loc_name must be specified")
    }
    if (!all(sapply(clones$data, function(clone) loc_name %in% colnames(clone@data)))){
      stop(paste0(loc_name, " must be a column in clones' data"))
    }
  }

  # --- Downsample Each Clone --- #
  clones$data <- lapply(clones$data, function(clone) downsampleSpaceClone(clone, downsample_map, method, x, y, loc_name, clone@clone, ...))
  clones$seqs <- sapply(clones$data, function(clone)nrow(clone@data))
  return(clones)
}


#' Downsample a clone in space according to a downsample map
#' 
#' @param    clone              airrClone object
#' @param    downsample_map     Data frame mapping all points onto their respective downsampled centroids, the output of \link{buildDownsampleMap}
#' @param    method             Either "coords" to match by x and y coordinates, or "locs" to match by location name (default "coords")
#' @param    x                  If method is "coords", the name of the column in \code{f_clone} corresponding to x coordinates
#' @param    y                  If method is "coords", the name of the column in \code{f_clone} corresponding to y coordinates
#' @param    loc_name               If method is "locs", the name of the column in \code{f_clone} corresponding to location names
#' @param    remove_areas       If specified, a data frame of areas to remove (created with \link{defineRemoveAreas}) to exclude from downsampling clone
#' 
#' @return   An \code{airrClone} object with locations downsampled according to the downsample map
#' 
#' @export
#'
downsampleSpaceClone <- function(clone, 
                              downsample_map,
                              method,
                              x = NULL,
                              y = NULL,
                              loc_name = NULL,
                              clone_id = NULL,
                              remove_areas = NULL) {

  # --- Check Parameters --- #
  if (!inherits(clone, "airrClone")){
    stop("clone must be an airrClone object")
  }

  if (!is.data.frame(downsample_map) | !(all(c("loc.x", "loc.y", "loc_name", "downsampled_centroid_name", "centroid.x", "centroid.y", "centroid_loc_name", "area", "key_area_name") %in% colnames(downsample_map)))){
    stop("downsample_map must be a data frame with columns loc.x, loc.y, loc_name, downsampled_centroid_name, centroid.x, centroid.y, centroid_loc_name, area, and key_area_name" )
  }

  if (!method %in% c("coords", "locs")){
    stop("method must be either 'coords' or 'locs'")
  }

  if (method == "coords"){
    if (is.null(x) | is.null(y)){
      stop("If method is 'coords', x and y must be specified")
    }
    if (!all(sapply(f_clones$data, function(clone) all(c(x, y) %in% colnames(clone@data))))){
      stop(paste0(x, " and ", y, " must be columns in f_clones' data"))
    }
  } else if (method == "locs"){
    if (is.null(loc_name)){
      stop("If method is 'locs', loc_name must be specified")
    }
    if (!all(sapply(f_clones$data, function(clone) loc_name %in% colnames(clone@data)))){
      stop(paste0(loc_name, " must be a column in f_clones' data"))
    }
  }

  if (!is.null(remove_areas)){
    if (!is.data.frame(remove_areas) | !(all(c("loc_name", "x", "y") %in% colnames(remove_areas)))){
      stop("If input, remove_areas must be a data frame with columns loc_name, x, y, and polygon_id")
    }
    remove_areas$polygon_id <- paste0("remove_area_", remove_areas$polygon_id)
  }

  if(is.null(clone_id)){clone_id <- ""}

  # --- Downsample the Clone --- #
  clone_data <- clone@data

  if (method == "coords"){
    # Harmonize types for joining
    clone_data[[x]] <- as.numeric(clone_data[[x]])
    clone_data[[y]] <- as.numeric(clone_data[[y]])

    # Extract locations of points in downsample_map
    downsample_locs <- downsample_map %>% dplyr::select(loc.x, loc.y)
    # Include remove_areas if provided
    if (!is.null(remove_areas)){
      remove_locs <- remove_areas %>% dplyr::select(x, y) %>% dplyr::rename(loc.x = x, loc.y = y)
      downsample_locs <- downsample_locs %>% dplyr::bind_rows(remove_locs)
    }

    # Match locations based on minimizing distance to spatial locations in downsample_map
    # Initialize vector to hold index of nearest location in downsample_map
    nearest_idx <- integer(nrow(clone_data))
    
    # Loop through locations in clone_data to find nearest location in downsample_map
    for (i in 1:nrow(clone_data)) {
      dx <- downsample_locs$loc.x - clone_data[[x]][i]
      dy <- downsample_locs$loc.y - clone_data[[y]][i]
      dists <- dx^2 + dy^2       # squared Euclidean distance
      nearest_idx[i] <- which.min(dists)
    
      # If remove_areas is specified, don't include any sequences that map to those areas
      if (!is.null(remove_areas)){
        if (nearest_idx[i] > nrow(downsample_map)) {
            nearest_idx[i] <- NA
        }
      }
    }

    # Join by nearest indices
    downsample_map_to_merge <- downsample_map[nearest_idx, , drop = FALSE]
    clone_data <- cbind(clone_data, downsample_map_to_merge)

    # Remove any sequences with NA nearest indices (i.e., mapped to remove_areas)
    numNAs <- sum(is.na(nearest_idx))
    if (numNAs > 0){
        warning(paste0(numNAs, "/", nrow(clone_data), " sequences in clone ", clone_id, " mapped to remove_areas and will be removed from the clone"))
    }
    clone_data <- clone_data[!is.na(nearest_idx), ]

  } else if (method == "locs"){
    
    # Harmonize nomenclature
    clone_data$loc_name <- as.character(clone_data[[loc_name]])
    clone_data <- dplyr::left_join(clone_data, downsample_map, by = "loc_name")

    # Remove any sequences that don't map to a location in downsample_map
    numNAs <- sum(is.na(clone_data$downsampled_centroid_name))
    if (numNAs > 0){
      warning(paste0(numNAs, "/", nrow(clone_data), " sequence locs are not in downsample_map and will be removed from the clone"))
      if (numNAs == nrow(clone_data)){
        stop(paste0("No locations in clone ", clone_id, " column ", loc_name, " match to downsample_map$loc_name. Cannot downsample clone, please remove and try again."))
      }
    }

    # Remove any sequences that don't map to a location in downsample_map
    clone_data <- clone_data %>% dplyr::filter(!is.na(downsampled_centroid_name))
  }

  # Update clone data and return
  clone@data <- clone_data
  return(clone)
}


#' Downsample the spatial coordinates while preserving key areas and removing remove areas
#' 
#' @param    cartesian_coords               Spatial cartesian coordinates data frame containing x and y columns and loc_names as row names
#' @param    key_areas                      key area data frame (created with \link{defineKeyAreas})
#' @param    remove_areas                   remove area data frame (created with \link{defineRemoveAreas}), or NULL if no areas to remove
#' @param    downsample_key_area            Downsample key areas too? (default TRUE)
#' @param    num_downsample_key_area        Number of centroids to downsample key area points onto (default 100)
#' @param    num_downsample_non_key         Number of centroids to downsample non-key area points onto (default 100)
#' @param    num_neighbors_non_boundary     Number of neighboring points needed to be considered non-boundary - depends on spatial sequencing method (default 6 for 10X Visium)
#' @param    neighbor_distance_threshold    Distance threshold to consider points as neighbors - calculated unless specified
#' @param    save_downsample_plot_path      Path to save downsampling plot, if NULL (default) then plot is not saved. Saves to PDF.
#' @param    seed_key_area                  Random seed for downsampling key area points
#' @param    seed_non_key                   Random seed for downsampling non-key area points
#'
#' @return   Data frame containing all points mapped to their downsampled centroids with columns loc.x, loc.y, loc_name, downsampled_centroid_name, centroid.x, centroid.y, centroid_loc_name, area (downsampled, centroid, or key_area)
#' 
#' @export
downsampleSpace <- function(cartesian_coords,
                              key_areas,
                              remove_areas = NULL,
                              downsample_key_area = TRUE,
                              num_downsample_key_area = 100,
                              num_downsample_non_key = 100,
                              num_neighbors_non_boundary = 6,
                              neighbor_distance_threshold = NULL,
                              save_downsample_plot_path = "downsample_plot.pdf",
                              seed_key_area = 123456,
                              seed_non_key = 7891011,
                              ...){
    
  # --- Check Parameters --- #
  if (!is.data.frame(cartesian_coords) | !(all(c("x", "y") %in% colnames(cartesian_coords)))){
    stop("cartesian_coords must be a data frame with columns x and y")
  }
  cartesian_coords$loc_name <- rownames(cartesian_coords)
  cartesian_coords <- cartesian_coords %>% dplyr::mutate(dplyr::across(c(x, y), ~ round(.x, 3))) # round coordinates to avoid floating point issues

  if (!is.data.frame(key_areas) | !(all(c("loc_name", "x", "y", "polygon_id") %in% colnames(key_areas)))){
    stop("key_areas must be a data frame with columns loc_name, x, y, polygon_id")
  }
  key_areas$polygon_id <- paste0("key_area_", key_areas$polygon_id)
  key_areas <- key_areas %>% dplyr::mutate(dplyr::across(c(x, y), ~ round(.x, 3))) # round coordinates to avoid floating point issues
  # Ensure that the key areas are a subset of the cartesian coords
  if (!all(do.call(paste0, key_areas %>% dplyr::select(x, y, loc_name)) %in% do.call(paste0, cartesian_coords))){
      stop("key_areas must be a subset of the coordinates in cartesian_coords")
  }
  
  if (!is.null(remove_areas)){
    if (!is.data.frame(remove_areas) | !(all(c("loc_name", "x", "y") %in% colnames(remove_areas)))){
      stop("If input, remove_areas must be a data frame with columns loc_name, x, y, and polygon_id")
    }
    remove_areas$polygon_id <- paste0("remove_area_", remove_areas$polygon_id)
    # Ensure that the remove areas are a subset of the cartesian coords
    if (!all(do.call(paste0, remove_areas %>% dplyr::select(x, y, loc_name)) %in% do.call(paste0, cartesian_coords))){
        stop("remove_areas must be a subset of the coordinates in cartesian_coords")
    }
  }

  if (downsample_key_area){
    if (!is.numeric(num_downsample_key_area) | num_downsample_key_area <= 0 | num_downsample_key_area > nrow(key_areas) | num_downsample_key_area %% 1 != 0){
      stop("num_downsample_key_area must be an integer greater than 0 and less than or equal to the number of points in key_areas")
    }
  }

  if (!is.numeric(num_downsample_non_key) | num_downsample_non_key <= 0 | num_downsample_non_key > (nrow(cartesian_coords) - nrow(key_areas)) | num_downsample_non_key %% 1 != 0){
    stop("num_downsample_non_key must be an integer greater than 0 and less than or equal to the number of non-key area points")
  }

  if (!is.numeric(num_neighbors_non_boundary) | num_neighbors_non_boundary <= 0 | num_neighbors_non_boundary %% 1 != 0){
    stop("num_neighbors_non_boundary must be a positive integer")
  }

  if (!is.null(neighbor_distance_threshold)){
    if (!is.numeric(neighbor_distance_threshold) | neighbor_distance_threshold <= 0){
      stop("If input, neighbor_distance_threshold must be a positive number")
    }
  }

  if (!is.null(save_downsample_plot_path)){
    if (!is.character(save_downsample_plot_path) | !endsWith(save_downsample_plot_path, ".pdf")){
      stop("If input, save_downsample_plot_path must be a character string representing the file path to save the plot and must end with .pdf")
    }
  }

  # --- Remove areas if requested --- #
  if (!is.null(remove_areas)){
    cartesian_coords <- cartesian_coords %>%
      dplyr::anti_join(remove_areas %>% dplyr::select(loc_name), by = "loc_name")
  }

  # --- Downsample the non-key area points --- #
  # Calculate the neighbor distance threshold if not provided
  if (is.null(neighbor_distance_threshold)){
    neighbor_threshold <- calculate_neighbor_distance_threshold(cartesian_coords %>% dplyr::anti_join(key_areas %>% dplyr::select(loc_name), by = "loc_name"),
                                                                 num_neighbors_non_boundary)
    if (is.null(neighbor_threshold)){
      stop("Could not calculate neighbor distance threshold automatically, please provide a value for neighbor_distance_threshold")
    }
  } else {
    neighbor_threshold <- neighbor_distance_threshold
  }
    
  # Get the centroids for the non key area
  centroids <- createCentroidsFPS(spatial_grid = cartesian_coords,
                                    exclude_points = key_areas %>% dplyr::select(x, y, loc_name),
                                    num_centroids = num_downsample_non_key,
                                    num_neighbors_non_boundary = num_neighbors_non_boundary,
                                    neighbor_distance_threshold = neighbor_threshold,
                                    seed = seed_non_key)
    
  # Build the downsampling map for the non key area (note that it includes all points)
  downsample_map <- buildDownsampleMap(spatial_grid = cartesian_coords,
                                           exclude_points = key_areas %>% dplyr::select(x, y, loc_name),
                                           centroids = centroids,
                                           rw_padding = neighbor_threshold)

  # --- Downsample the key area points if requested --- #
  if (downsample_key_area){
    # Get the centroids for the key area
    key_area_centroids <- createCentroidsFPS(spatial_grid = key_areas %>% dplyr::select(x, y, loc_name),
                                               num_centroids = num_downsample_key_area,
                                               num_neighbors_non_boundary = num_nearby_for_non_boundary_point,
                                               neighbor_distance_threshold = distance_for_nearby,
                                               seed = seed_key_area)
      
    # Build the downsampling map for the key area
    downsample_map_key_area <- buildDownsampleMap(spatial_grid = key_areas %>% dplyr::select(x, y, loc_name),
                                                    centroids = key_area_centroids,
                                                    rw_padding = neighbor_threshold)
      
    downsample_map_key_area <- downsample_map_key_area %>%
                                   dplyr::mutate(downsampled_centroid_name = paste0("key_area_", downsampled_centroid_name)) %>%
                                   dplyr::mutate(area = recode(area, "centroid" = "key_area_centroid"))
          
    # Replace the key areas in the original map with this one
    downsample_map <- downsample_map %>% dplyr::filter(area != "excluded_from_downsampling") %>%
        dplyr::bind_rows(downsample_map_key_area)
    downsample_map$area <- factor(downsample_map$area, levels = c("key_area_centroid", "centroid", "downsampled"))

    
    # Plot the downsampling if requested
    if (!is.null(save_downsample_plot_path)){
      all_centroids <- key_area_centroids %>%
          dplyr::mutate(centroid_name = paste0("key_area_", centroid_name)) %>%
          dplyr::bind_rows(centroids)

      g <- plot_downsample(spatial_grid = cartesian_coords,
                             centroids = all_centroids,
                             downsample_map = downsample_map)
      # TODO: Add outline to key areas maybe?  That would be nice.
    
      pdf(save_downsample_plot_path, width = 8, height = 6)
        print(g)
      dev.off()
    }
  }
  
  # --- Prepare Outputs --- #
  # Add the key area names to the downsampled mapping, NA if not in a key area
  downsample_map <- downsample_map %>% mutate(key_area_name = key_areas$polygon_id[match(loc_name, key_areas$loc_name)])

  if (!downsample_key_area){
    # Fix the area names for the key areas if we didn't downsample key areas
    downsample_map <- downsample_map %>%
                         dplyr::mutate(downsampled_centroid_name = gsub("exclude_point", "key_area_point_", downsampled_centroid_name),
                                         area = recode(area, "excluded_from_downsampling" = "key_area"))

    # Plot the downsampling if requested
    if (!is.null(save_downsample_plot_path)){
      g <- plot_downsample(spatial_grid = cartesian_coords,
                             key_areas = key_areas,
                             centroids = centroids,
                             downsample_map = downsample_map)
    
      pdf(save_downsample_plot_path, width = 8, height = 6)
        print(g)
      dev.off()
    }  
  }

  # --- Compute adjacency and contextual matrices if we decide to go that route --- #  
  # TODO: Implement later if needed
    # Get the downsampled adjacency matrix
    # adjacency_matrix <- compute_adjacency_matrix(spatial_grid = new_coords,
    #                                              mapping = downsample_map, 
    #                                              rw_padding = rw_padding, 
    #                                              plot_adjacency = FALSE)
    
    # # Get the downsampled key area contextual matrix
    # contextual_matrix <- construct_contextual_matrix2(mapping = downsample_map,
    #                                                  adjacency_matrix = adjacency_matrix)
      
    # out_list <- list("downsampled_mapping" = downsample_map,
    #                 "adjacency_matrix" = adjacency_matrix,
    #                 "contextual_matrix" = contextual_matrix)
    
    # return(out_list)

  return(downsample_map)
}



######################################################################################################
#######                                                                                         ######
#######                                    Helper Functions                                     ######
#######                                                                                         ######
######################################################################################################

#' Determine a distance threshold via kernel density estimation to consider points as neighbors
#' 
#' @param    cartesian_coords               Spatial cartesian coordinates data frame containing x and y columns and loc_names as row names
#' @param    num_neighbors_non_boundary     Number of nearby points needed to be considered non-boundary 
#' 
#' @return  Distance threshold to consider points as neighbors
#' 
calculate_neighbor_distance_threshold <- function(cartesian_coords, num_neighbors_non_boundary){
  # Calculate distance matrix
  dist_matrix <- stats::dist(cartesian_coords %>% dplyr::select(x, y), method = "euclidean") %>% as.matrix()

  # For each point, get the distance to its num_neighbors_non_boundary-th nearest neighbor
  neighbor_distances <- apply(dist_matrix, 1, function(row) {
    sorted_distances <- sort(row[row > 0])  # Exclude distance to self
    return(sorted_distances[num_neighbors_non_boundary])
  })
  
  # Remove any NAs, NaNs, or Infs
  neighbor_distances <- neighbor_distances[!is.na(neighbor_distances) & !is.nan(neighbor_distances) & !is.infinite(neighbor_distances)]

  # Need at least 3 unique distances for density estimation
  unique_distances <- unique(neighbor_distances)
  if (length(unique_distances) < 3 ) {
    stop("The automated distance threshold calculation requires at least 3 unique distance values.\n",
         "Found distance values: ", paste(unique_distances, collapse=", "))
  }
  
  # Estimate ideal bandwidth
  bandwidth <- kedd::h.ucv(unique_distances, 4)$h
  
  # Kernel density estimation
  dens <- KernSmooth::bkde(neighbor_distances, bandwidth = bandwidth, canonical = TRUE)
  xdens <- dens$x
  ydens <- dens$y
  
  # Find distance threshold
  tryCatch(threshold <- xdens[which(diff(sign(diff(ydens))) == 2)[1] + 1], 
           error = function(e) {
             warning('No minimum was found between two modes.')
             return(NULL) })
  
  return(threshold)
}

#' Utility function to count the number of decimal places in a number
#' @param    x         Vector of numeric values
#'
#' @return   Number of decimal places in each value of x
#' 
count_decimals <- function(x) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }

  dec_vec <- sapply(x, function(val) {
    if ((val %% 1) != 0) {
      nchar(strsplit(sub('0+$', '', as.character(val)), ".", fixed=TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  })

  return(dec_vec)
}

#' Define polygon(s) from a spatial grid by encapsulating all points that lie on or within the previously defined boundary
#' 
#' @param    spatial_grid         Data frame with columns x, y, and loc_name of the spatial grid
#' @param    by_all_or_boundary   Either "all" or "boundary" to specify if all or boundary points of the polygon are provided
#' @param    by_names_or_coords   Either "names" or "coords" to specify if the polygon is defined by location names or coordinates
#' @param    names_or_coords      Vector of location names or data frame with columns x and y defining the polygon
#' 
#' @return  Subset of spatial_grid data frame that comprises the points in the polygon
#' 
define_polygon <- function(spatial_grid, by_all_or_boundary, by_names_or_coords, names_or_coords){

  # --- Check Parameters --- #
  if (!is.data.frame(spatial_grid) | !(all(colnames(spatial_grid) %in% c("x", "y", "loc_name")))){
    stop("spatial_grid must be a data frame with columns x, y, and loc_name")
  }
  spatial_grid <- spatial_grid %>% dplyr::select(x, y, loc_name) %>%
    dplyr::mutate(x = round(x, 3), y = round(y, 3))

  if (!(by_all_or_boundary == "all" | by_all_or_boundary == "boundary")){
    stop("by_all_or_boundary must be either 'all' or 'boundary'")
  }

  if (!(by_names_or_coords == "names" | by_names_or_coords == "coords")){
    stop("by_names_or_coords must be either 'names' or 'coords'")
  }

  if (by_names_or_coords == "names"){
    if (!is.vector(names_or_coords)){
      stop("names_or_coords must be a vector of names")
    }
    if (!all(names_or_coords %in% spatial_grid$loc_name)){
      stop("all names must be in spatial_grid")
    }
  } else {
    if (!is.data.frame(names_or_coords) | !(all(colnames(names_or_coords) %in% c("x", "y")))){
      stop("names_or_coords must be a data frame with columns x and y")
    }
    if (!all(names_or_coords$x %in% spatial_grid$x) | !all(names_or_coords$y %in% spatial_grid$y)){
      stop("all coordinates must be in spatial_grid")
    }
  }

  # --- Define the polygon when all points are specified --- #
  if (by_all_or_boundary == "all"){
    if (by_names_or_coords == "names"){
      key_area <- spatial_grid %>% filter(loc_name %in% names_or_coords)
    } else {
      key_area <- spatial_grid %>% filter(x %in% names_or_coords$x & y %in% names_or_coords$y)
    }
  }

  # --- Define the polygon when just the boundary is specified --- #
  else if (by_all_or_boundary == "boundary"){
    if (by_names_or_coords == "names"){

      # Make sure the last name is the same as the first name
      if (names_or_coords[1] != names_or_coords[length(names_or_coords)]){
        stop("The first and last names must be the same to define a boundary")
      }

      # Convert names to coordinates
      boundary_coords <- spatial_grid[match(names_or_coords, spatial_grid$loc_name), ] %>% select(x, y)

    } else if (by_names_or_coords == "coords"){

      # Make sure the last coordinate is the same as the first coordinate
      if (!all(names_or_coords[1, ] == names_or_coords[nrow(names_or_coords), ])){
        stop("The first and last coordinates must be the same to define a boundary")
      }

      boundary_coords <- names_or_coords
    }
    
    # Convert boundary points to an sf polygon
    boundary_sf <- sf::st_sfc(sf::st_polygon(list(as.matrix(boundary_coords))))

    # Convert spatial_grid to sf points
    spatial_sf <- sf::st_as_sf(spatial_grid, coords = c("x", "y"))

    # Identify points inside the polygon or on the boundary
    inside <- sf::st_within(spatial_sf, boundary_sf, sparse = FALSE) %>% as.vector()
    on_boundary <- sf::st_touches(spatial_sf, boundary_sf, sparse = FALSE) %>% as.vector()
    inside_or_on_boundary <- inside | on_boundary

    # Subset spatial_grid to only include inside points
    key_area <- spatial_grid[inside_or_on_boundary, ]
    
  }
  return(key_area)
}



#' Utility function to find the name of the column with a 1 in a given row of a matrix with 1s and 0s
#' @param    mat         Matrix with 1s and 0s
#' @param    rowname     Row name to search
#'
#' @return   Name of the column with a 1 in the specified row, or NA if none
#' 
find_one <- function(mat, rowname){
  mat_row <- mat[rowname, ]
  out <- names(which(mat_row == 1))
  # If there are no 1s, return an NA
  if (length(out) == 0){
    return(NA)
  } else {
    return(out)
  }
}



#' Plot gene expression in Cartesian space for defineKeyAreas
#' 
#' @param    seurat_obj         Spatial Seurat object containing coordinates and gene expression data
#' @param    gene               Gene to plot, must be in Features(seurat_obj)
#' @param    assay              Assay containing gene expression data (default Spatial)
#' @param    dataset            Dataset within the assay containing gene expression data (default scale.data)
#' @param    point_size         Size of points in the plot (default 3)
#'
#' @return   list with plotly object and data frame of expression data, for use in defineKeyAreas
#' 
plot_gene_expression_cartesian <- function(seurat_obj, gene, assay = "Spatial", dataset = "scale.data", point_size = 3, ...){

  # --- Check Parameters --- #
  if (!inherits(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }

  if (!assay %in% names(seurat_obj@assays)) {
    stop(paste0("Assay '", assay, "' not found in seurat_obj@assays."))
  }
  
  if (!dataset %in% slotNames(seurat_obj[[assay]])) {
    stop(paste0("Dataset '", dataset, "' not found in assay '", assay, "'."))
  }
    
  if (!gene %in% Features(seurat_obj)) {
    stop("gene must be a feature in seurat_obj")
  }

  # --- Get gene expression and spatial coordinates --- #
  # Get the expression of the gene
  gex <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = dataset)
  gene_expression <- data.frame(cell = names(gex[gene, ]), expression = gex[gene, ])

  # Get the coordinates of the spots
  coords <- Seurat::GetTissueCoordinates(seurat_obj)
  cartesian_coords <- convertCoordinates(coords, x="imagecol", y="imagerow", direction="p2c")
  cartesian_coords$cell <- rownames(cartesian_coords)

  gex_coords <- dplyr::left_join(cartesian_coords, gene_expression, by = c("cell" = "cell")) %>% dplyr::select(cell, x, y, expression)
  colnames(gex_coords) <- c("ID", "x", "y", "expression")

  # --- Create the plot --- #
  # Set plot dimensions to maintain aspect ratio
  width = 800
  height = width * (max(gex_coords$y) - min(gex_coords$y)) / (max(gex_coords$x) - min(gex_coords$x))

  # Plot the expression
  g <- ggplot2::ggplot() + 
          ggplot2::geom_point(
            data = gex_coords,
            ggplot2::aes(x = .data$x, y = .data$y, color = .data$expression),
            size = point_size
          ) +
          ggplot2::theme_bw() +
          ggplot2::theme(panel.grid = element_blank()) +
          ggplot2::scale_color_viridis_c(...)
  gply <- plotly::ggplotly(g, source = "plot_click") %>% plotly::layout(autosize = F, width = width, height = height)
  return(list(plot = gply, data = gex_coords))
}

#' Plot downsampling for downsampleSpace
#' 
#' @param    spatial_grid        Data frame with columns x, y, and loc_name covering the spatial grid
#' @param    key_areas           Data frame with columns x, y, and loc_name covering the key areas (if not downsampled as well)
#' @param    centroids           Data frame with columns x, y, loc_name, and centroid_name covering the centroids
#' @param    downsample_map      Data frame with columns loc.x, loc.y, loc_name, downsampled_centroid_name, centroid.x, centroid.y, centroid_loc_name, and area (from buildDownsampleMap)
#' @param    point_size          Size of points in the plot (default 2)
#'
#' @return   list with plotly object and data frame of expression data, for use in defineKeyAreas
#' 
plot_downsample <- function(spatial_grid, key_areas = NULL, centroids, downsample_map, point_size = 2){

  # --- Check Parameters --- #
  if (!is.data.frame(spatial_grid) | !(all(c("x", "y", "loc_name") %in% colnames(spatial_grid)))){
    stop("spatial_grid must be a data frame with columns x, y, and loc_name")
  }
  spatial_grid <- spatial_grid %>% select(x, y, loc_name)
  
  if (!is.null(key_areas)){
    if (!is.data.frame(key_areas) | !(all(c("x", "y", "loc_name") %in% colnames(key_areas)))){
      stop("key_areas must be a data frame with columns x, y, and loc_name")
    }
    key_areas <- key_areas %>% select(x, y, loc_name)
    # Ensure key_areas are a subset of spatial_grid
    if (!all(do.call(paste0, key_areas) %in% do.call(paste0, spatial_grid))){
      stop("key_areas must be a subset of spatial_grid")
    }
  }
    
  if (!is.data.frame(centroids) | !(all(c("x", "y", "loc_name", "centroid_name") %in% colnames(centroids)))){
    stop("centroids must be a data frame with columns x, y, loc_name, and centroid_name")
  }
  centroids <- centroids %>% select(x, y, loc_name, centroid_name)
  # Ensure that centroids are a subset of spatial_grid
  if (!all(do.call(paste0, centroids[, c("x", "y", "loc_name")]) %in% do.call(paste0, spatial_grid))){
    stop("centroids must all be in spatial_grid")
  }

  if (!is.data.frame(downsample_map) | !(all(c("loc.x", "loc.y", "loc_name", "downsampled_centroid_name", "centroid.x", "centroid.y", "centroid_loc_name", "area") %in% colnames(downsample_map)))){
    stop("downsample_map must be a data frame with columns loc.x, loc.y, loc_name, downsampled_centroid_name, centroid.x, centroid.y, centroid_loc_name, and area.\n,
         Please use the output of buildDownsampleMap.")
  }

  # Make a column in spatial_grid to indicate what type of area a point lives in and another column for the centroid name (if a point is a key area or itself a centroid then loc+name == centroid_loc_name)
  spatial_grid$area <- "base"
  if (!is.null(key_areas)){
    spatial_grid$area[spatial_grid$loc_name %in% key_areas$loc_name] <- "key_area"
  }
  spatial_grid$area[spatial_grid$loc_name %in% centroids$loc_name] <- "centroid"
  spatial_grid$area <- factor(spatial_grid$area, levels = c("base", "centroid", "key_area"))
  spatial_grid$loc_name <- as.character(spatial_grid$loc_name)
  spatial_grid <- dplyr::left_join(spatial_grid, downsample_map[, c("loc_name", "centroid_loc_name")], by = c("loc_name"))
  spatial_grid$centroid_loc_name <- factor(spatial_grid$centroid_loc_name)

  # Set up a contrasting color palette to improve visualization
  n_centroids <- spatial_grid %>% filter(area == "centroid") %>% nrow()
  quotient <- n_centroids %/% 32
  remainder <- n_centroids %% 32
  palette <- c(rep(pals::glasbey(32), quotient), pals::glasbey(remainder))
  names(palette) <- spatial_grid %>% filter(area == "centroid") %>% pull(centroid_loc_name)


  # Plot the spatial grid with highlighted key areas and centroids
  g <- ggplot2::ggplot() + 
         ggplot2::geom_point(data = filter(spatial_grid, area == "centroid"), aes(x=x, y=y, fill = centroid_loc_name), color = "black", size = point_size, shape = 21) +
         ggplot2::geom_point(data = filter(spatial_grid, area == "base"), aes(x=x, y=y, color = centroid_loc_name), size = point_size/2) +
         ggplot2::scale_color_manual(values = palette, guide = "none") +
         ggplot2::scale_fill_manual(values = palette, guide = "none") +
         ggplot2::theme_bw() +
         ggplot2::theme(panel.grid = element_blank(), legend.position = "none") 
  if (!is.null(key_areas)){
      g <- g + ggplot2::geom_point(data = filter(spatial_grid, area == "key_area"), aes(x=x,y=y), size = point_size, fill = "green3", color = "black", shape = 21) 
  }
  return(g)
}


#' Set up the UI for the defineKeyAreas Shiny app
#'
#' @return   UI for the Shiny app
#' 
set_up_ui <- function(){
  ui <- shiny::fluidPage(
    shiny::titlePanel("Define Polygons by Drawing Boundaries"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::textInput("gene", "Gene Name", value = "AICDA"),
        shiny::numericInput("point_size", "Point Size", value = 2, min = 1, max = 10),
        shiny::textInput("save_name", "Filename", value = "polygons.csv"),
        shiny::actionButton("remove_last_point", "Remove Last Point"),
        shiny::actionButton("reset_points", "Reset Selection"),
        shiny::actionButton("new_polygon", "New Polygon"),
        shiny::actionButton("save_points", "Save and Exit"),
        shiny::tableOutput("selected_points"),
        width = 5
      ),
      shiny::mainPanel(plotly::plotlyOutput("gene_plot", width = "100%"), width = 7)
    )
  )
  return(ui)
}