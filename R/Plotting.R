#Plotting functions for trees

# Blend set of colors
# 
# \code{combineColors} blends colors
# @param    x    vector of states
# @param    pal  The colorbrewer palette to use
#
# @return   A color hex code representing the average of input colors
#
combineColors <- function(x,pal){
    cols <- rowMeans(grDevices::col2rgb(pal[as.character(x)]))
    col <- grDevices::rgb(cols[1],cols[2],cols[3],maxColorValue=255)
    return(col)
}

#' Get a color palette for a predefined set of trait values
#' 
#' \code{getPalette} Gets a color palette for a predefined set of trait values
#' @param    states   states in model
#' @param    palette  The colorbrewer palette to use
#'
#' @return   A named vector with each state corresponding to a color
#'
#' @seealso \link{getTrees}, \link{plotTrees}
#' @export
getPalette <- function(states,palette){
    if(palette == "AmG"){
                #M        D         G13       A1        G24
        pal <- c("#000000", "#696969", "#33a02c", "#1f78b4", "#e31a1c",
            #E         A2        G         A
            "#dd3497", "#6a3d9a", "#33a02c", "#1f78b4")
        names(pal) <- c("M", "D", "G31", "A1", "G24", "E", "A2", "G", "A") 
    }else if(palette == "FullIg"){
                #M        D         G3        G1        A1        G2
        pal <- c("#000000", "#696969", "#b15928", "#33a02c", "#1f78b4", "#e31a1c",
            #G4       #E         #A2
            "#ff7f00", "#dd3497", "#6a3d9a")
        names(pal) <- c("IgM", "IgD", "IgG3", "IgG1", "IgA1", "IgG2", "IgG4",
            "IgE", "IgA2") 
    }else if (palette == "IgAmG" || palette == "IgAmGA"){
                #M        D         G13       A1        G24
        pal <- c("#000000", "#696969", "#33a02c", "#1f78b4", "#e31a1c",
            #E         A2        G         A
            "#dd3497", "#6a3d9a", "#33a02c", "#1f78b4")
        names(pal) <- c("IgM", "IgD", "IgG31", "IgA1", "IgG24", "IgE", "IgA2", 
            "IgG", "IgA") 
    }else{
        pal <- RColorBrewer::brewer.pal(n=length(unique(states)), name=palette)
        names(pal) <- as.character(unique(states))
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
condenseTrees <- function(trees, states, palette){
    if(class(trees) == "phylo"){
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
    cv <- unlist(lapply(combs,function(x)combineColors(x,palette)))
    margl <- unlist(lapply(combs,function(x)paste(x,collapse=",")))
    nt$node.label <- margl[(tipn+1):noden]
    nt$node.color <- cv
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
        palette <- c(palette,"ambig"="#808080")
    }
    for(n in 1:length(trees)){
        nt <- trees[[n]]
        combs <- strsplit(nt$state, split=",")
        if(ambig == "blend"){
            cv <- unlist(lapply(combs, function(x)combineColors(x, palette)))
        }else if(ambig == "grey"){
            combs[unlist(lapply(combs,function(x)length(x) > 1))] <- "ambig"
            cv <- unlist(lapply(combs, function(x)combineColors(x, palette)))
            nt$state <- unlist(combs)
        }else{
            stop("ambig parameter not specified")
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
#' @param    trees        A tibble containing \code{phylo} and \code{airrClone}
#'                        objects
#' @param    nodes        color internal nodes if possible?
#' @param    tips         color tips if possible?
#' @param    tipsize      size of tip shape objects
#' @param    scale        width of branch length scale bar
#' @param    node_palette color palette for nodes
#' @param    tip_palette  color palette for tips
#' @param    layout       rectangular or circular tree layout?
#' @param    node_nums    plot internal node numbers?
#' @param    tip_nums     plot tip numbers?
#' @param    title        use clone id as title?
#' @param    labelsize    text size
#' @param    base         recursion base case (don't edit)
#'
#' @return   a grob containing a tree plotted by \code{ggtree}.
#'
#' @details
#' Function uses \code{ggtree} functions to plot tree topologlies estimated by 
#' \link{getTrees}, and \link{bootstrapTrees}. Object can be further modified with 
#' \code{ggtree} functions. Please check out 
#' https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html and
#' cite \code{ggtree} in addition to \code{dowser} if you use this function.
#'  
#' @seealso \link{getTrees}, \link{bootstrapTrees}
#' @examples
#' data(ExampleClones)
#' trees <- getTrees(ExampleClones[10,])
#' plotTrees(trees)[[1]]
#' @export
plotTrees <- function(trees, nodes=FALSE, tips=NULL, tipsize=NULL, 
    scale=0.01, node_palette="Dark2", tip_palette=node_palette, base=FALSE,
    layout="rectangular", node_nums=FALSE, tip_nums=FALSE, title=TRUE,
    labelsize=NULL){

    tiptype = "character"
    if(!base){
        cols <- c()
        if(!is.null(tips) && nodes && sum(tip_palette != node_palette) == 0){
            tipstates <- unique(c(unlist(lapply(trees$data,function(x)
                unique(x@data[[tips]])))))
            if(is.numeric(tipstates)){
                stop("Can't currently plot numeric tip values and node values")
            }
            tipstates = c(tipstates,"Germline")
            nodestates <- unique(unlist(lapply(trees$trees,function(x)
                    unique(unlist(strsplit(x$state,split=",")))
                    )))
            combpalette <- getPalette(c(nodestates,tipstates),node_palette)
            trees$trees <- colorTrees(trees$trees,palette=combpalette)
            nodestates <- unlist(lapply(trees$trees,function(x){
                colors <- x$node.color
                names(colors) <- x$state
                colors
                }))
            nodepalette <- nodestates[unique(names(nodestates))]
            cols <- c(combpalette,nodepalette[!names(nodepalette) %in% names(combpalette)])
        }else{
            if(!is.null(tips)){
                tipstates <- unique(c(unlist(lapply(trees$data,function(x)
                    unique(x@data[[tips]])))))
                if(is.numeric(tipstates)){
                    tiptype <- "numeric"
                    cols <- range(tipstates)
                }else{
                    tipstates = c(tipstates,"Germline")
                    if(is.null(names(tip_palette))){
                        tip_palette <- getPalette(tipstates,tip_palette)
                        tip_palette <- tip_palette[!is.na(names(tip_palette))]
                    }else{
                        nfound <- tipstates[!tipstates %in% names(tip_palette)]
                        if(length(nfound) > 0){
                            stop(paste(nfound,"not found in tip_palette"))
                        }
                    }
                    cols <- tip_palette
                }
            }
            if(nodes){
                if(is.null(names(node_palette))){
                    nodestates <- unique(unlist(lapply(trees$trees,function(x)
                        unique(unlist(strsplit(x$state,split=",")))
                        )))
                    statepalette <- getPalette(nodestates,node_palette)
                    statepalette <- statepalette[!is.na(names(statepalette))]
                }else{
                    statepalette <- node_palette
                }
                trees$trees <- colorTrees(trees$trees,palette=statepalette)
                
                nodestates <- unlist(lapply(trees$trees,function(x){
                    colors <- x$node.color
                    names(colors) <- x$state
                    colors
                    }))
                nodepalette <- nodestates[unique(names(nodestates))]
                cols <- c(tip_palette,nodestates)
            }
        }
        
        ps <- lapply(1:nrow(trees),function(x)plotTrees(trees[x,],
            nodes=nodes,tips=tips,tipsize=tipsize,scale=scale,node_palette=node_palette,
            tip_palette=tip_palette,base=TRUE,layout=layout,node_nums=node_nums,
            tip_nums=tip_nums,title=title,labelsize=labelsize))
        if(!is.null(tips) || nodes){
            ps  <- lapply(ps,function(x){
                    x <- x + theme(legend.position="right",
                    legend.box.margin=margin(0, -10, 0, 0))+
                    guides(color=guide_legend(title="State"))
                    if(tiptype == "character"){
                        x <- x + scale_color_manual(values=cols)
                    }else{
                        x <- x + scale_color_distiller(limits=cols,
                            palette=tip_palette)
                    }})
        }
        return(ps)
    }

    tree <- trees$trees[[1]]
    data <- trees$data[[1]]
    p <- ggtree::ggtree(tree,layout=layout)
    if(!is.null(data)){
        if(class(data) != "list"){
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
            if(class(data@data[[n]]) == "numeric" || class(data@data[[n]]) == "integer"){
                gl[[n]] <- NA
            }else if(class(data@data[[n]]) == "character"){
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
    if(!is.null(tips)){
        if(is.null(data)){
            stop("dataframe must be provided when tip trait specified")
        }
        if(!is.null(tipsize)){
            if(class(tipsize) == "numeric"){
                p <- p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips)),size=tipsize)
            }else if(class(tipsize) == "character"){
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
#' @param  nrow    number of rows per page
#' @param  ncol    size of tip shape objects
#' @param  ...    optional arguments passed to grDevices::pdf
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

