#Plotting functions for parsimony-labeled trees

#blend set of colors
combineColors = function(x,pal){
	cols = rowMeans(grDevices::col2rgb(pal[x]))
	col = grDevices::rgb(cols[1],cols[2],cols[3],maxColorValue=255)
	return(col)
}

#' Get a color palette for a predefined set of trait values
#' 
#' \code{getPalette} Gets a color palette for a predefined set of trait values
#' @param    palette  The colorbrewer palette to use
#' @param    states    states in model
#'
#' @return   A named vector with each state corresponding to a color
#'
#' @seealso \link{getTrees}, \link{plotTrees}
#' @export
getPalette = function(palette,states){
	if(palette == "AmG"){
				#M        D         G13       A1        G24
		pal = c("#000000","#696969","#33a02c","#1f78b4","#e31a1c",
			#E         A2        G         A
			"#dd3497","#6a3d9a","#33a02c","#1f78b4")
		names(pal) = c("M","D","G31","A1","G24","E","A2","G","A") 
	}else if(palette == "FullIg"){
				#M        D         G3        G1        A1        G2
		pal = c("#000000","#696969","#b15928","#33a02c","#1f78b4","#e31a1c",
			#G4       #E         #A2
			"#ff7f00","#dd3497","#6a3d9a")
		names(pal) = c("IgM","IgD","IgG3","IgG1","IgA1","IgG2","IgG4",
			"IgE","IgA2") 
	}else if (palette == "IgAmG" || palette == "IgAmGA"){
				#M        D         G13       A1        G24
		pal = c("#000000","#696969","#33a02c","#1f78b4","#e31a1c",
			#E         A2        G         A
			"#dd3497","#6a3d9a","#33a02c","#1f78b4")
		names(pal) = c("IgM","IgD","IgG31","IgA1","IgG24","IgE","IgA2","IgG","IgA") 
	}else{
		pal = RColorBrewer::brewer.pal(n = length(unique(states)), name = palette)
		names(pal) = unique(states)
	}
	return(pal)
}

#' Condense a set of equally parsimonious node labels into a single tree
#' 
#' \code{condenseTrees} Condenses a set of equally parsimonious node labels into a single tree
#' @param    trees     List of the same tree with equally parsimonious labels
#' @param    palette   Named vector with a color per state
#' @param    states    States in model
#'
#' @return   a \code{phylo} object representing all represented internal node states
#'
#' @export
condenseTrees = function(trees,states,palette){
	if(class(trees) == "phylo"){
		trees = list(trees)
		class(trees) = "multiPhylo"
	}
	if(is.null(names(states))){
		names(states) = states
	}
	nt = trees[[1]]
	tipn = length(trees[[1]]$tip.label)
	noden = 2*tipn-1
	combs = list()
	for(i in 1:(noden)){
		combs[[i]] = sort(unique(states[
			unlist(lapply(trees,function(x){
				lab=x$node.comment[i]
				if(!lab %in% names(states)){
					stop(paste(lab,"not found in names of state vector"))
				}
				lab
				}))]))
	}
	nt$unique = unique(unlist(combs))
	cv = unlist(lapply(combs,function(x)combineColors(x,palette)))
	margl = unlist(lapply(combs,function(x)paste(x,collapse=",")))
	nt$node.label = margl[(tipn+1):noden]
	nt$node.color = cv
	nt$state = margl
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
#'
#' @seealso \link{getPalette}, \link{getTrees}, \link{plotTrees}
#' @export
colorTrees <- function(trees,palette,ambig="blend"){
    ntrees <- list()
    if(ambig == "grey"){
    	palette = c(palette,"ambig"="#808080")
    }
    for(n in 1:length(trees)){
    	nt <- trees[[n]]
		tipn <- length(nt$tip.label)
		noden <- 2*tipn-1
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
		ntrees[[n]] = nt
    }
    #class(ntrees) <- "multiPhylo"
    return(ntrees)
}

#' Plot a tree with colored internal node labels using ggtree
#' 
#' \code{plotTrees} plots a tree or group of trees
#' @param    trees      A tibble containing \code{phylo} and \code{changeoClone} objects
#' @param    data     	(optional) 
#' @param    nodes   	color internal nodes if possible?
#' @param    tips 		color tips if possible?
#' @param    trait    	trait to use to color the tips
#' @param 	 tipsize 	size of tip shape objects
#' @param 	 data 		list of \code{changeoClone} objects used to generate \code{tree}
#' @param  	 scale 		width of branch length scale bar
#'
#' @return   a grob containing a tree plotted by \code{ggtree}.
#'
#' @details
#' Function uses \code{ggtree} functions to plot tree topologlies estimated by \link{getTrees},
#' and \link{bootstrapTrees}. Object can be further modified with \code{ggtree} functions. Please
#' check out https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html and
#' cite \code{ggtree} in addition to \code{dowser} if you use this function.
#'  
#' @seealso \link{getTrees}, \link{bootstrapTrees}
#' @examples
#' \dontrun{
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb,trait="sample_id")
#'
#' trees = getTrees(clones[1:2])
#' plotTrees(trees[[1]])
#' }
#' @export
#' # make tips the trait arugment
#' # make separate nodes and tips color scale
#' # make tip size variable
plotTrees = function(trees,data=NULL,nodes=TRUE,tips=NULL,tipsize=NULL,scale=0.01,
	node_palette="Dark2",tip_palette=node_palette,base=FALSE,layout="rectangular"){
	if(!base){
		if(!is.null(tips) && nodes && tip_palette == node_palette){
			tipstates = unique(unlist(lapply(trees$DATA,function(x)unique(x@data[[tips]]))))
			nodestates = unique(unlist(lapply(trees$TREE,function(x)
					unique(unlist(strsplit(x$state,split=",")))
					)))
			combpalette = getPalette(node_palette,c(nodestates,tipstates))
			trees$TREE = colorTrees(trees$TREE,palette=combpalette)
			nodestates = unlist(lapply(trees$TREE,function(x){
				colors = x$node.color
				names(colors) = x$state
				colors
				}))
			nodepalette = nodestates[unique(names(nodestates))]
			cols = c(combpalette,nodepalette[!names(nodepalette) %in% names(combpalette)])
		}else{
			if(!is.null(tips)){
				tipstates = unique(unlist(lapply(trees$DATA,function(x)unique(x@data[[tips]]))))
				if(is.atomic(tip_palette)){
					tip_palette = getPalette(tip_palette,tipstates)
					tip_palette = tip_palette[!is.na(names(tip_palette))]
				}else{
					nfound = tipstates[!tipstates %in% names(tip_palette)]
					if(length(nfound) > 0){
						stop(paste(nfound,"not found in tip_palette"))
					}
				}
			}
			if(nodes){
				nodestates = unique(unlist(lapply(trees$TREE,function(x)
					unique(unlist(strsplit(x$state,split=",")))
					)))
				statepalette = getPalette(node_palette,nodestates)
				statepalette = statepalette[!is.na(names(statepalette))]
				trees$TREE = colorTrees(trees$TREE,palette=statepalette)
				
				nodestates = unlist(lapply(trees$TREE,function(x){
					colors = x$node.color
					names(colors) = x$state
					colors
					}))
				nodepalette = nodestates[unique(names(nodestates))]
			}
			cols = c(tip_palette,nodestates)
		}
		
		ps = lapply(1:nrow(trees),function(x)plotTrees(trees[x,],
			nodes=nodes,tips=tips,tipsize=tipsize,scale=scale,node_palette=node_palette,
			tip_palette=tip_palette,base=TRUE,layout=layout))
		
		ps  = lapply(ps,function(x)
				x = x + theme(legend.position="right",
		    	legend.box.margin=margin(0, -10, 0, 0))+
		    	scale_color_manual(values=cols)+
		    	guides(color=guide_legend(title="State")))
		return(ps)
	}

	tree = trees$TREE[[1]]
	data = trees$DATA[[1]]
	p = ggtree::ggtree(tree,layout=layout)
	if(!is.null(data)){
		if(class(data) != "list"){
			data = list(data)
		}
		index = which(unlist(lapply(data,function(x)x@clone == tree$name)))
		if(length(index) == 0){
			stop("clone",tree$name," not found in list of clone objects")
		}
		if(length(index) > 1){
			stop("clone",tree$name," found more than once in list of clone objects")
		}
		data = data[[index]]
		p = p %<+% data@data
	}
	if(!is.null(tree$pars_recon)){
		if(nodes){
			p = p + aes(color=tree$state)
		}
	}
	if(!is.null(tips)){
		if(is.null(data)){
			stop("dataframe must be provided when tip trait specified")
		}
		if(!is.null(tipsize)){
			if(class(tipsize) == "numeric"){
				p = p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips)),size=tipsize)
			}else if(class(tipsize) == "character"){
				p = p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips),
					size=!!rlang::sym(tipsize)))
			}
		}else{
			p = p + ggtree::geom_tippoint(aes(color=!!rlang::sym(tips)))
		}
	}
	p = p + ggtree::geom_treescale(width=scale)
	p
}

# get position grid to arrange ggtree plots in a 
# grid
arrangeGrid = function(trees,transform="sqrt",
	pack=FALSE,widths=NULL,wratio=1){
	plots = list()
	if(is.null(widths)){
		widths = unlist(lapply(trees,
			function(x)length(x$tip.label)))
	}
	if(transform=="sqrt"){
		widths= sqrt(widths)
	}else if(transform=="cubic"){
		widths = widths^(1/3)
	}else{
		widths= widths
	}
	
	l = length(widths)
	mint = min(widths)
	block = ceiling(widths/mint)
	index = 1:l
	area = sum(block^2)
	minwidth = ceiling(sqrt(area/wratio))
	grid = matrix(NA,nrow=l,ncol=minwidth)
	
	for(i in 1:l){
		for(j in 1:minwidth){
			if(length(block) == 0){
				break
			}
			if(is.na(grid[i,j])){
				nal = sum(is.na(grid[i,j:minwidth]))
				if(sum(!is.na(grid[i,j:minwidth])) > 0){
					stopper = which(!is.na(grid[i,j:minwidth]))[1]
					stopper = stopper + j -1
					nal = sum(is.na(grid[i,j:stopper]))
				}
				if(nal < block[1] && !pack){
					next
				}
				opt = block - nal
				opt[opt>0] = -Inf
				block_index = which.max(opt)
				icoord = i + block[block_index]-1
				jcoord = j + block[block_index]-1
				for(it in i:icoord){
					for(jt in j:jcoord){
						grid[it,jt] = index[block_index]
					}
				}
				block = block[-block_index]
				index = index[-block_index]
			}
		}
	}
	grid = grid[rowSums(is.na(grid)) != minwidth,]

	if(sum(sqrt(table(grid))-floor(sqrt(table(grid)))) != 0){
		print("grid construction failed :-(")
		print(sqrt(table(grid)))
		stop()
	}
	return(grid)
}
