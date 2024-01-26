#Plotting functions for trees

# Blend set of colors
# 
# \code{combineColors} blends colors
# @param    x    vector of states
# @param    pal  The colorbrewer palette to use
#
# @return   A color hex code representing the average of input colors
#
combineColors <- function(x, pal){
    x <- as.character(x)
    if(sum(is.na(pal[x])) != 0){
        stop(paste(x[!x %in% names(pal)]," not included in palette!"))
    }
    cols <- rowMeans(grDevices::col2rgb(pal[as.character(x)]))
    col <- grDevices::rgb(cols[1],cols[2],cols[3],maxColorValue=255)
    return(col)
}

#' Get a color palette for a predefined set of trait values. 
#' 'Germline' defaults to black unless specified.
#' 
#' \code{getPalette} Gets a color palette for a predefined set of trait values
#' @param    states   states in model
#' @param    palette  The colorbrewer palette to use
#'
#' @return   A named vector with each state corresponding to a color
#'
#' @seealso \link{getTrees}, \link{plotTrees}
#' @export
getPalette <- function(states, palette){
  # CGJ 12/6/23
  if(length(palette) == 1 && palette == "AmG"){
    #M        D         G13       A1        G24
    pal <- c("#000000", "#696969", "#33a02c", "#1f78b4", "#e31a1c",
             #E         A2        G         A
             "#dd3497", "#6a3d9a", "#33a02c", "#1f78b4")
    names(pal) <- c("M", "D", "G31", "A1", "G24", "E", "A2", "G", "A") 
  }else if(length(palette) == 1 && palette == "FullIg"){
    #M        D         G3        G1        A1        G2
    pal <- c("#000000", "#696969", "#b15928", "#33a02c", "#1f78b4", "#e31a1c",
             #G4       #E         #A2
             "#ff7f00", "#dd3497", "#6a3d9a")
    names(pal) <- c("IgM", "IgD", "IgG3", "IgG1", "IgA1", "IgG2", "IgG4",
                    "IgE", "IgA2") 
  }else if (length(palette) == 1 && palette == "IgAmG" || length(palette) == 1 &&
            palette == "IgAmGA"){
    #M        D         G13       A1        G24
    pal <- c("#000000", "#696969", "#33a02c", "#1f78b4", "#e31a1c",
             #E         A2        G         A
             "#dd3497", "#6a3d9a", "#33a02c", "#1f78b4")
    names(pal) <- c("IgM", "IgD", "IgG31", "IgA1", "IgG24", "IgE", "IgA2", 
                    "IgG", "IgA") 
  }else{
    # 12/20/23 CGJ -- changed the if elses to not freak out over the named vector
    # also updated this section to check for a named vector or just a RBrewer input
    if(length(palette) > 1 && !is.null(names(palette))){
      # check that all states (besides Germline) are found in palette
        pal_names <- names(palette)
        pal <- palette
        if(sum(!states %in% pal_names) != 0){
            disjoint <- states[!states %in% pal_names]
            if(length(disjoint) == 1 && disjoint == "Germline"){
                warning("Germline not included in palette, setting to black")
                pal["Germline"] <- "#000000"
            }else{
                stop(paste("States not in palette, please add:",
                    paste(disjoint[disjoint != "Germline"], collapse=", ")))
            }
        }
    } else{
      # changed 11/14/23 CGJ
      # test if the palette is too small for what they are trying to do
      nongerm_states <- unique(states[states != "Germline"])
      nstates <- dplyr::n_distinct(nongerm_states)
      pal_test <- suppressWarnings(tryCatch(
        RColorBrewer::brewer.pal(n=nstates, name=palette),
        error=function(e)e))
      if(nstates > length(pal_test)){
        # if it is send a warning and replace the palette
        warning(paste("There are", nstates, "unique tips specified",
                      "which is more than the", palette, "allows. Switching to a",
                      "larger palette."))
        if(nstates > 69){
          stop(paste("There are", nstates, "unique states in a specified tip",
                     "plotting variable. There are more states than what can be plotted."))
        }
        # this finds all the quantitative colors in RColorBrewer
        qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
        # get the palette of all of them together
        pal <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        pal <- unique(pal)
        # remove the bright yellow right at the beginning
        pal <- pal[!pal %in% '#FFFF99']
        
        # cut to where you need
        pal <- pal[1:nstates]
        names(pal) <- as.character(nongerm_states)
        if("Germline" %in% states){
            pal["Germline"] <- "#000000"
        }
      } else{
        pal <- RColorBrewer::brewer.pal(n=nstates, name=palette)
        names(pal) <- as.character(nongerm_states)
        if("Germline" %in% states){
            pal["Germline"] <- "#000000"
        }
      }
    }
  }
  return(pal)
}

#' Condense a set of equally parsimonious node labels into a single tree
#' 
#' \code{condenseTrees} Condenses a set of equally parsimonious node labels 
#' into a single tree
#' @param    trees     List of the same tree with equally parsimonious labels
#' @param    palette   Named vector with a color per state
#' @param    states    States in model
#'
#' @return   a \code{phylo} object representing all represented internal node states
#'
#' @export
condenseTrees <- function(trees, states, palette=NULL){
    if(!is.null(palette)){
        warning("palette option is deprecated in condenseTrees, specify in plotTrees")
    }
    if(is(trees,"phylo")){
        trees <- list(trees)
        class(trees) <- "multiPhylo"
    }
    if(is.null(names(states))){
        names(states) <- states
    }
    nt <- trees[[1]]
    tipn <- length(trees[[1]]$tip.label)
    noden <- 2*tipn-1
    combs <- list()
    for(i in 1:(noden)){
        combs[[i]] <- sort(unique(states[
            unlist(lapply(trees,function(x){
                lab=x$node.comment[i]
                if(!lab %in% names(states)){
                    stop(paste(lab,"not found in names of state vector"))
                }
                lab
                }))]))
    }
    nt$unique <- unique(unlist(combs))
    #cv <- unlist(lapply(combs,function(x)combineColors(x,palette)))
    margl <- unlist(lapply(combs,function(x)paste(x,collapse=",")))
    nt$node.label <- margl[(tipn+1):noden]
    #nt$node.color <- cv
    nt$state <- margl
    return(nt)
}

#' Get a color palette for a predefined set of trait values
#' 
#' \code{colorTree} Gets a color palette for a predefined set of trait values
#' @param    trees   list of phylo objects with assigned internal node states
#' @param    palette named vector of colors (see \link{getPalette})
#' @param    ambig   how should ambiguous states be colored (blend or grey)
#'
#' @return   A list of colored trees
#' 
#' @details Trees must have node states represented in a "states" vector. By default,
#' ambiguous states (separated by ",") have their colors blended. If 
#' 
#' @seealso \link{getPalette}, \link{getTrees}, \link{plotTrees}
#' @export
colorTrees <- function(trees, palette, ambig="blend"){
    ntrees <- list()
    if(ambig == "grey"){
        palette <- c(palette,"Ambiguous"="#808080")
    }
    for(n in 1:length(trees)){
        nt <- trees[[n]]
        combs <- strsplit(nt$state, split=",")
        if(ambig == "blend"){
            cv <- unlist(lapply(combs, function(x)combineColors(x, palette)))
        }else if(ambig == "grey"){
            combs[unlist(lapply(combs,function(x)length(x) > 1))] <- "Ambiguous"
            cv <- unlist(lapply(combs, function(x)combineColors(x, palette)))
            nt$state <- unlist(combs)
        }else{
            stop("ambig parameter must be either 'blend' or 'grey'")
        }
        nt$node.color <- cv
        ntrees[[n]] <- nt
    }
    #class(ntrees) <- "multiPhylo"
    return(ntrees)
}

#' Plot a tree with colored internal node labels using ggtree
#' 
#' \code{plotTrees} plots a tree or group of trees
#' @param    trees              A tibble containing \code{phylo} and \code{airrClone}
#'                              objects
#' @param    nodes              color internal nodes if possible?
#' @param    tips               color tips if possible?
#' @param    tipsize            size of tip shape objects
#' @param    scale              width of branch length scale bar
#' @param    palette            color palette for tips and/or nodes. Can supply a named vector
#'                              for all tip states, or a palette named passed to
#'                              ggplot2::scale_color_brewer (e.g. "Dark2", "Paired", "Set1") or
#'                              ggplot2::scale_color_distiller (e.g. RdYlBu) or
#' @param    common_scale       strecth plots so branches are on same scale?
#'                              determined by sequence with highest divergence
#' @param    layout             rectangular or circular tree layout?
#' @param    node_nums          plot internal node numbers?
#' @param    tip_nums           plot tip numbers?
#' @param    title              use clone id as title?
#' @param    labelsize          text size
#' @param    base               recursion base case (don't edit)
#' @param    ambig              How to color ambiguous node reconstructions? (grey or blend)
#' @param    bootstrap_scores    Show bootstrap scores for internal nodes? See getBootstraps.
#' @param    height_intervals   plot 95 percent HPD intervals (only BEAST trees)
#' @param    densitree          plot posterior density trees (only BEAST trees)
#' @param    densitree_sample   trees to plot in densitree
#' @param    node_palette       deprecated, use palette
#' @param    tip_palette        deprecated, use palette
#' @param    height_width       width of height intervals 
#' @param    height_color       color of height intervals 
#'
#' @return   a grob containing a tree plotted by \code{ggtree}.
#'
#' @details
#' Function uses \code{ggtree} functions to plot tree topologlies estimated by 
#' \link{getTrees}, and \link{findSwitches}. Object can be further modified with 
#' \code{ggtree} functions. Please check out 
#' https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html and
#' cite \code{ggtree} in addition to \code{dowser} if you use this function.
#'  
#' @seealso \link{getTrees}, \link{findSwitches}
#' @examples
#' data(ExampleClones)
#' trees <- getTrees(ExampleClones[10,])
#' plotTrees(trees)[[1]]
#' @export
plotTrees <- function(trees, nodes=FALSE, tips=NULL, tipsize=NULL, 
    scale=0.01, palette="Dark2", base=FALSE,
    layout="rectangular", node_nums=FALSE, tip_nums=FALSE, title=TRUE,
    labelsize=NULL, common_scale=FALSE, ambig="grey", bootstrap_scores=FALSE,
    height_intervals=FALSE, tip_palette=NULL, node_palette=NULL, height_width=NULL,
    height_color="red", densitree=FALSE, densitree_sample=100){

    tiptype = "character"
    # CGJ 12/12/23 add check to see if the color palettes are unnamed vectors 
    if(is.null(names(palette)) && length(palette) > 1){
      stop("palette must be either a named vector or ColorBrewer palette name (string)")
    }
    # KBH 12/21/23 Deprecate tip_palette and node_palette, now all is palette
    if(!base && (!is.null(tip_palette) || !is.null(node_palette))){
        warning("tip_palette and node_palette are deprecated, please use palette instead.")
        if(!is.null(tip_palette)){
            palette <- tip_palette
        }else{
            palette <- node_palette
        }
    }
    if(height_intervals && densitree){
        stop("Cannot plot both height_intervals and densitree")
    }
    if(!base){
        cols <- c()
        # set up global tip and node palette
        if(!is.null(tips) && nodes){
            tipstates <- unique(c(unlist(lapply(trees$data,function(x)
                unique(x@data[[tips]])))))
            if(is.numeric(tipstates)){
                stop("Can't currently plot numeric tip values and node values")
            }
            tipstates = c(sort(tipstates),"Germline")
            if(is.null(names(palette))){
              nodestates <- sort(unique(unlist(lapply(trees$trees,function(x)
                unique(unlist(strsplit(x$state,split=",")))
              ))))
            } else{
              nodestates <- names(palette)
            }
            combpalette <- getPalette(unique(c(nodestates,tipstates)),palette)
            trees$trees <- colorTrees(trees$trees,palette=combpalette,ambig=ambig)
            nodestates <- unlist(lapply(trees$trees,function(x){
                colors <- x$node.color
                names(colors) <- x$state
                colors
                }))
            nodepalette <- nodestates[unique(names(nodestates))]
            cols <- c(combpalette,nodepalette[!names(nodepalette) %in% names(combpalette)])
        }else if(!is.null(tips)){
            # set up global tip palette
            if(!is.null(tips)){
                tipstates <- unique(c(unlist(lapply(trees$data,function(x)
                    unique(x@data[[tips]])))))
                if(is.numeric(tipstates)){
                    tiptype <- "numeric"
                    cols <- range(tipstates)
                }else{
                    tipstates = c(sort(tipstates),"Germline")
                    if(is.null(names(palette))){
                        palette <- getPalette(tipstates,palette)
                        palette <- palette[!is.na(names(palette))]
                    }else{
                        palette <- getPalette(tipstates,palette)
                        nfound <- tipstates[!tipstates %in% names(palette)]
                        if(length(nfound) > 0){
                            stop(paste(nfound,"not found in palette"))
                        }
                    }
                    cols <- palette
                }
            }
        }else if(nodes){
            # set up global node palette
            if(is.null(names(palette))){
                nodestates <- unique(unlist(lapply(trees$trees,function(x)
                    unique(unlist(strsplit(x$state,split=",")))
                    )))
                statepalette <- getPalette(sort(nodestates),palette)
                statepalette <- statepalette[!is.na(names(statepalette))]
            }else{
                statepalette <- palette
            }
            trees$trees <- colorTrees(trees$trees,palette=statepalette, ambig=ambig)
            
            nodestates <- unlist(lapply(trees$trees,function(x){
                colors <- x$node.color
                names(colors) <- x$state
                colors
                }))
            nodepalette <- nodestates[unique(names(nodestates))]
            cols <- c(palette,nodepalette[!names(nodepalette) %in% names(palette)])
        }
        if(common_scale){
            # get maximum divergence value
            max_div <- max(unlist(lapply(trees$trees, function(x)max(getDivergence(x)))))
        }
        
        ps <- lapply(1:nrow(trees),function(x)plotTrees(trees[x,],
            nodes=nodes,tips=tips,tipsize=tipsize,scale=scale,palette=palette, node_palette=node_palette,
            tip_palette=tip_palette,base=TRUE,layout=layout,node_nums=node_nums,
            tip_nums=tip_nums,title=title,labelsize=labelsize, ambig=ambig, 
            bootstrap_scores=bootstrap_scores, height_intervals=height_intervals,
            height_width=height_width, height_color=height_color, densitree=densitree,
            densitree_sample=densitree_sample))
        if(!is.null(tips) || nodes){
            ps  <- lapply(ps,function(x){
                    x <- x + theme(legend.position="right",
                    legend.box.margin=margin(0, -10, 0, 0))+
                    guides(color=guide_legend(title="State"))
                    if(tiptype == "character"){
                        x <- x + scale_color_manual(values=cols)
                    }else{
                        x <- x + scale_color_distiller(limits=cols,
                            palette=palette)
                    }})
        }
        if(common_scale){
             ps  <- lapply(ps,function(x){
                x <- x + xlim(0, max_div*1.05)
            })
        }
        return(ps)
    }

    tree <- trees$trees[[1]]
    data <- trees$data[[1]]
    p <- ggtree::ggtree(tree,layout=layout)

    if(densitree){
        phylos <- tree$tree_posterior
        if(length(phylos) == 0){
            stop("Tree posterior not found. Be sure tree_posterior=TRUE in getTrees and/or readBEAST")
        }
        p <- p + ggtree::geom_tiplab()
        pdata <- p$data
        pdata <- pdata[order(pdata$y),]
        tip_order <- pdata[pdata$isTip,]$label

        phylos <- phylos[floor(seq(1, length(phylos), length=densitree_sample))]
        p <- ggtree::ggdensitree(phylos, tip.order=tip_order, 
            alpha=2/densitree_sample)
    }

    #add bootstrap scores to ggplot object
    if(bootstrap_scores){
        scores <- unlist(lapply(tree$nodes, function(x)x$bootstrap_value))
        if(is.null(scores)){
            stop(paste("No bootstrap scores found in tree",tree$name))
        }
        p$data$bootstrap_score = scores[p$data$node]
    }
    if(!is.null(data)){
        if(!is(data,"list")){
            data <- list(data)
        }
        index <- which(unlist(lapply(data,function(x)x@clone == tree$name)))
        if(length(index) == 0){
            stop("clone",tree$name," not found in list of clone objects")
        }
        if(length(index) > 1){
            stop("clone",tree$name," found more than once in list of clone objects")
        }
        data <- data[[index]]
        gl <- dplyr::tibble(sequence_id="Germline")
        for(n in names(data@data)){
            if(is(data@data[[n]],"numeric") || is(data@data[[n]], "integer")){
                gl[[n]] <- NA
            }else if(is(data@data[[n]], "character")){
                gl[[n]] <- "Germline"
            }else{
                gl[[n]] <- "Germline"
            }
        }
        data@data <- rbind(data@data,gl)
        p <- p %<+% data@data
    }
    if(!is.null(tree$pars_recon)){
        if(nodes){
            p <- p + aes(color=tree$state)
        }
    }
    if(!is.null(tree$tree_method)){
        if(grepl("BEAST", tree$tree_method) && !densitree){
            intervals <- lapply(tree$nodes, function(x)x$height_0.95_HPD)
            node_data <- dplyr::tibble(node=1:length(intervals),
              height_0.95_HPD=intervals)
            p <- p %<+% node_data

            if(height_intervals){
                if(length(intervals) == 0){
                    stop("height_0.95_HPD not found in nodes, was tree built with BEAST?")
                }
                if(!is.null(height_width)){
                    p <- p + ggtree::geom_range("height_0.95_HPD", color=height_color,
                        size=height_width, alpha=.5) 
                }else{
                    p <- p  + ggtree::geom_range("height_0.95_HPD", color=height_color, 
                        alpha=.5) 
                }
            }
        }
    }
    if(!is.null(tips)){
        if(is.null(data)){
            stop("dataframe must be provided when tip trait specified")
        }
        if(!is.null(tipsize)){
            if(is(tipsize, "numeric")){
                p <- p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips)),size=tipsize)
            }else if(is(tipsize, "character")){
                p <- p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips),
                    size=!!rlang::sym(tipsize)))
            }
        }else{
            p <- p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips)))
        }
    }
    if(scale != FALSE){
        p <- p + ggtree::geom_treescale(width=scale)
    }
    if(title){
        p <- p + ggtitle(data@clone)
    }
    if(bootstrap_scores){
        if(is.null(labelsize)){
            p <- p + ggtree::geom_label(data=p$data[!p$data$isTip,],
                aes(label=!!rlang::sym("bootstrap_score")),label.padding = unit(0.1, "lines"),
                label.size=0.1)
        }else{
            p <- p + ggtree::geom_label(data=p$data[!p$data$isTip,],
                aes(label=!!rlang::sym("bootstrap_score")),label.padding = unit(0.1, "lines"),
                label.size=0.1,size=labelsize)
        }
    }
    if(node_nums){
        if(is.null(labelsize)){
            p <- p + ggtree::geom_label(data=p$data[!p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1)
        }else{
            p <- p + ggtree::geom_label(data=p$data[!p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1,size=labelsize)
        }
    }
    if(tip_nums){
        if(is.null(labelsize)){
            p <- p + ggtree::geom_label(data=p$data[p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1)
        }else{
            p <- p + ggtree::geom_label(data=p$data[p$data$isTip,],
                aes(label=!!rlang::sym("node")),label.padding = unit(0.1, "lines"),
                label.size=0.1,size=labelsize)
        }
    }
    p
}

#' Simple function for plotting a lot of trees into a pdf
#' 
#' \code{treesToPDF} exports trees to a pdf in an orderly fashion
#' @param  plots list of tree plots (from plotTrees)
#' @param  file  output file name
#' @param  nrow  number of rows per page
#' @param  ncol  number of columns per page
#' @param  ...   optional arguments passed to grDevices::pdf
#' 
#' @return   a PDF of tree plots
#'  
#' @seealso \link{plotTrees}
#' @examples
#' \dontrun{
#' data(ExampleClones)
#' trees <- getTrees(ExampleClones[10,])
#' plots <- plotTrees(trees)
#' treesToPDF(plots,"test.pdf",width=5,height=6)
#' }
#' @export
treesToPDF = function(plots, file, nrow=2, ncol=2, ...){
    treepage = nrow*ncol
    rm = treepage - length(plots) %% treepage
    if(rm != treepage){
        for(i in (length(plots)+1):(length(plots)+rm)){
            plots[[i]] = ggplot(data.frame())
        }
    }
    s = seq(1,length(plots),by=treepage)
    grDevices::pdf(file,...)
    for(start in s){
        gridExtra::grid.arrange(grobs=
            plots[start:(start + treepage -1)], ncol=ncol)    
    }
    grDevices::dev.off()
}

#' Simple function for plotting BEAST trace parameter traces
#' 
#' \code{plotTraces} exports trees to a pdf in an orderly fashion
#' @param  clones output from getTrees using BEAST
#' @param  burnin percent of initial samples to discard (1-100)
#' @param  file   file name for printing plots
#' @param  width  width of plot in inches
#' @param  height height of plot in inches
#' @param  ess    add ESS to facets?
#' @param  ...   optional arguments passed to grDevices::pdf
#' 
#' @return   a PDF of trace plots
#'  
#' @seealso \link{getTrees}
#' @export
plotTraces = function(clones, burnin=10, file=NULL, width=8.5, height=11, ess=TRUE, ...){
    if(burnin > 100 || burnin < 0){
        stop("burnin must be between 0 and 100")
    }
    plots <- list()
    for(i in 1:nrow(clones)){
        post <- clones$trees[[i]]$parameters_posterior
        if(nrow(post) == 0){
            stop(paste("parameters_posterior empty for clone", clones$clone_id[[i]]))
        }
        sample_range = range(post$Sample)
        burn <- floor(burnin/100 * (sample_range[2] - sample_range[1]))
        post <- filter(post, !!rlang::sym("Sample") >= burn)

        if(ess){
            essv = post %>%
                group_by(parameter) %>%
                summarize(ess = mcmcse::ess(value))

            labels = paste0(essv$parameter, " | ESS: ", floor(essv$ess))
            names(labels) = essv$parameter
        }else{
            labels = unique(post$parameter)
            names(labels) = labels
        }

        post$parameter <- factor(post$parameter, levels=unique(post$parameter))
        plots[[i]] <-  ggplot(post, aes(x=Sample,y=value)) + geom_line() + 
            facet_wrap(parameter~., scales="free_y", ncol=2,
                labeller = labeller(parameter=labels)) +
            theme_bw() + ggtitle(clones$clone_id[[i]]) + ylab("Parameter value")
    }

    if(!is.null(file)){
        grDevices::pdf(file, width=width, height=height)
        for(i in 1:length(plots)){
            gridExtra::grid.arrange(grobs=plots[i], ncol=1)    
        }
        grDevices::dev.off()
    }else{
        return(plots)
    }
}

#' Simple function for plotting Bayesian skyline plots
#' 
#' \code{plotSkylines} Simple Bayesian skyline plots
#' @param  clones output from getTrees using BEAST
#' @param  file   pdf file name for printing plots
#' @param  width  width of plot in inches if file specified
#' @param  height height of plot in inches if file specified
#' @param  ess    add ESS to facets?
#' @param  ...   optional arguments passed to grDevices::pdf
#' 
#' @return   if no file specified, a list of ggplot objects. If file specified
#' will plot to specified file
#'  
#' @seealso \link{getSkylines} \link{readBEAST} \link{getTrees}
#' @export
plotSkylines = function(clones, file=NULL, width=8.5, height=11){
    plots <- list()
    for(i in 1:nrow(clones)){
        skyline <- clones$skyline[[i]]
        if(is.na(clones$skyline[i])){
            warning(clones$clone_id[i], "skyline not found, skipping")
            next
        }

        plots[[i]] <- ggplot(skyline, aes(x=bin, y=median, ymin=lci, ymax=uci)) +
            geom_ribbon(fill = "grey70") + scale_y_log10() + theme_bw() + 
            geom_line() + xlab("Time") + ylab("Effective pop. size") +
            ggtitle(clones$clone_id[i])
    }

    if(!is.null(file)){
        grDevices::pdf(file, width=width, height=height)
        for(i in 1:length(plots)){
            gridExtra::grid.arrange(grobs=plots[i], ncol=1)    
        }
        grDevices::dev.off()
    }else{
        return(plots)
    }
}








