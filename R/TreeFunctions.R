# Functions for performing discrete trait analysis on B cell lineage trees

#Write a clone's sequence alignment to a fasta file
writeFasta = function(c,fastafile,germid,trait=NULL,dummy=FALSE){
	clone = c@clone
	append = FALSE
	if(!is.null(trait)){
		c@data$SEQUENCE_ID = paste(c@data$SEQUENCE_ID,c@data[,trait],sep="_")
	}
	for(i in 1:nrow(c@data)){
		if(i > 1){append=TRUE}
		write(paste0(">",c@data[i,]$SEQUENCE_ID),
			file=fastafile,append=append)
		if(!dummy){
			write(c@data[i,]$SEQUENCE,file=fastafile,append=TRUE)
		}else{
			write("ATG",file=fastafile,append=TRUE)
		}
	}
	write(paste0(">",germid),file=fastafile,append=append)
	if(!dummy){
		write(c@germline,file=fastafile,append=TRUE)
	}else{
		write("ATG",file=fastafile,append=TRUE)
	}
	return(fastafile)
}

#' Read in a parsimony model file
#' 
#' \code{readModelFile} Filler
#' @param    file  		parimony model file.
#' @param 	 useambig 	use ambiguous naming as specified in the file?
#'
#' @return   A named vector containing the states of the model
#'
#' @seealso \link{makeModelFile}, \link{bootstrapTrees}, \link{getTrees}
#'
#' @export
readModelFile = function(file,useambig=FALSE){
	#set up pallete
	mfile = readLines(file)
	mstart = which(mfile == "#STATES")
	mend = which(mfile == "")
	mend = min(mend[which(mfile == "") > mstart])
	states = mfile[(mstart+1):(mend-1)]
	names(states) = states
	
	if(useambig){
		astart = which(mfile == "#AMBIGUOUS")
		ambigs = mfile[(astart+1):length(mfile)]
		asplit = strsplit(ambigs,split=" ")
		ambig = unlist(lapply(asplit,function(x)x[1]))
		names(ambig) = unlist(lapply(asplit,function(x)x[2]))
		ambig = ambig[names(ambig) != "M"]
		nambig = states[!states %in% names(ambig)]
		names(nambig) = nambig
		ambig = c(ambig,nambig)
	}else{
		ambig = states
	}
	return(ambig)
}

#' Make a parsimony model file
#' 
#' \code{makeModelFile} Filler
#' @param    file  			model file name to write.
#' @param 	 states 		vector of states to include in model.
#' @param    constraints 	constraints to add to model.
#'
#' @return   Name of model file
#'
#' @details
#' Currently the only option for \code{constraints} is "irrev", which
#' forbids switches moving from left to right in the \code{states} vector.
#'  
#' @seealso \link{readModelFile}, \link{getTrees}, \link{bootstrapTrees}
#'
#' @export
makeModelFile = function(file,states,constraints=NULL){
	write("#BEGIN", file=file)
	write(length(states), file=file, append=TRUE)
	write("", file=file, append=TRUE)
	write("#STATES", file=file, append=TRUE)
	for(a in states){
		write(a, file=file, append=TRUE)
	}
	write("", file=file, append=TRUE)
	write("#CONSTRAINTS", file=file, append=TRUE)
	if(!is.null(constraints)){
		if(constraints=="irrev"){
			for(i in 1:(length(states)-1)){
				for(j in (i+1):length(states)){
					write(paste(states[j],states[i],"1000"), 
						file=file, append=TRUE)
				}
			}
		}else{
			stop("Contraint not recognized.")
		}
	}
	write("", file=file, append=TRUE)
	write("", file=file, append=TRUE)
	write("#AMBIGUOUS", file=file, append=TRUE)
	for(a in states){
		write(paste("GERM",a), file=file, append=TRUE)
	}
	return(file)
}

#read in switches data set
readSwitches = function(file){
	t = read.table(file,sep="\t",stringsAsFactors=FALSE)
	names(t) = c("REP","FROM","TO","SWITCHES")
	switches = t
	switches
}

#remove uniformative columns from data and germline
cleanAlignment = function(clone,seq="SEQUENCE"){
	if(seq=="HLSEQUENCE"){
    	g = strsplit(clone@hlgermline[1],split="")[[1]]
    }else{
    	g = strsplit(clone@germline[1],split="")[[1]]
    }
    #g = g[1:(length(g) - length(g)%%3)]
    sk = strsplit(clone@data[[seq]],split="")
    sites=seq(1,length(g)-3,by=3)
    ns = c()
    for(i in sites){
        l=lapply(sk,function(x) paste(x[i:(i+2)],collapse="")=="NNN")
        ns = c(ns,sum(unlist(l)),sum(unlist(l)),sum(unlist(l)))
    }
    informative = ns != length(sk)
    l=lapply(sk,function(x) x=paste(x[informative],collapse=""))
    gm=paste(g[informative],collapse="")
    if(.hasSlot(clone,"region")){
    	clone@region = clone@region[informative]	
    }
    if(seq=="HLSEQUENCE"){
    	clone@hlgermline=gm
    }else{
    	clone@germline=gm
    }
    clone@data[[seq]]=unlist(l)
    return(clone)
}

#uncollapse data based on trait column
uncollapse = function(clone, trait){
	data = clone@data
	newdata = data.frame()
	for(i in 1:nrow(data)){
		row = data[i,]
		types = strsplit(row[,trait],split=",")[[1]]
		for(type in types){
			newrow = row
			newrow[,trait] = type
			#newrow$DUPCOUNT = row$DUPCOUNT/length(types)
			newrow$SEQUENCE_ID = paste0(newrow$SEQUENCE_ID,"_",type)
			newdata = rbind(newdata,newrow)
		}
	}
	clone@data = newdata
	clone
}

#make bootstrap replicate clones
bootstrapClones  = function(clone, reps=100){
	sarray = strsplit(clone@data$SEQUENCE,split="")
	garray = strsplit(clone@germline,split="")[[1]]
	index = 1:stats::median(nchar(clone@data$SEQUENCE))
	bootstraps = list()
	for(i in 1:reps){
		clone_copy = clone
		sindex = sample(index,length(index),replace=TRUE)
		#print(paste(length(unique(sindex)),length(unique(index))))
		clone_copy@data$SEQUENCE = unlist(lapply(sarray,
			function(x)paste(x[sindex],collapse="")))
		clone_copy@germline = paste(garray[sindex],collapse="")
		bootstraps[[i]] = clone_copy
	}
	bootstraps
}

#reconstruct using maximum parsimony implemented in igphyml
reconIgPhyML = function(file, modelfile, cloneid, 
	igphyml="igphymlp",	mode="switches", type="recon",
	nproc=1, quiet=0, rm_files=FALSE, rm_dir=NULL, 
	states=NULL, palette=NULL, resolve=2, rseed=NULL){

	igphyml <- path.expand(igphyml)
	if(file.access(igphyml, mode=1) == -1) {
        stop("The file ", igphyml, " cannot be executed.")
    }
    if(!file.exists(file)) {
        stop("The repertoire file ", file, " cannot be found.")
    }
    if(!file.exists(modelfile)) {
        stop("The model file ", modelfile, " cannot be found.")
    }
	if(!mode %in% c("switches","trees")){
		stop(paste("mode must be either switches or trees"))
	}
	if(!type %in% c("recon","permute","permuteAll")){
		stop(paste("type must be either recon, permute, or permuteAll"))
	}

	recon = paste0(file,"_igphyml_parstats_",type,".txt")
	logfile = paste0(file,".log")
	log = paste(">>",logfile)
	permute = ""
	if(type == "permute"){
		permute = "--permute"
	}
	if(type == "permuteAll"){
		permute = "--permuteAll"
	}
	if(is.null(rseed)){
		rseed = ""
	}else{
		rseed = paste("--r_seed",rseed)
	}
	command = paste("--repfile",file,
		"--recon",modelfile,"--threads",nproc,"--polyresolve",resolve,
		"-m HLP -o n --motifs WRC_2:0 --hotness 0 --run_id",type,permute,
		rseed,log)
	params = list(igphyml,command,stdout=TRUE,stderr=TRUE)
	if(quiet > 2){
		print(paste(params,collapse=" "))
	}
	status = tryCatch(do.call(base::system2, params), error=function(e){
		print(paste("igphyml error, trying again: ",e));
		cat(paste(readLines(logfile),"\n"))
		return(e)
		}, warning=function(w){
		print(paste("igphyml warnings, trying again: ",w));
		cat(paste(readLines(logfile),"\n"))
		return(w)
		})
	if(length(status) != 0){
		status = tryCatch(do.call(base::system2, params), error=function(e){
		print(paste("igphyml error, again! quitting: ",e));
		cat(paste(readLines(logfile),"\n"))
		stop()
		}, warning=function(w){
		print(paste("igphyml warnings, again! quitting: ",w));
		cat(paste(readLines(logfile),"\n"))
		stop()
		})
	}

	if(quiet > 2){
		cat(paste(readLines(logfile),"\n"))
	}
	if(mode == "switches"){
		recons = readSwitches(recon)
		recons$CLONE = cloneid
		recons$TYPE = toupper(type)
		results = dplyr::as_tibble(recons)
	}else{
		if(is.null(states)){
			states = readModelFile(modelfile)
		}
		if(is.null(palette)){
			palette = getPalette("Dark2",states)
		}
		results = readLineages(file,states,palette,"recon",quiet)
		results = lapply(results,function(x){
			x$pars_recon="igphyml";
			x})
	}
	if(rm_files){
		lines = readLines(file)
		for(i in 2:length(lines)){
			temp = strsplit(lines[i],split="\t")[[1]]
			unlink(paste0(temp[1],"*"))
			unlink(paste0(temp[2],"*"))
		}
		unlink(paste0(file,"*"))
	}
	if(!is.null(rm_dir)){
		if(quiet > 1){
			print(paste("rming dir",rm_dir))
		}
		unlink(rm_dir,recursive=TRUE)
	}
	return(results)
}

#read in all trees from a nexus file
readLineages = function(file, states=NULL, palette="Dark2",
	run_id="", quiet=TRUE, append=NULL, format="nexus"){
	trees = list()
	switch = data.frame()
	t = readLines(file)
	if(length(t) == 0){
		return(list())
	}
	for(i in 2:length(t)){
		fasta = strsplit(t[i],split="\t")[[1]][1]
		if(quiet > 0){print(fasta)}
		if(is.null(append)){
			if(run_id == ""){
				tf = suppressWarnings(phylotate::read_annotated(
					paste0(fasta,"_igphyml_jointpars.nex"),
					format=format))
			}else{
				tf = suppressWarnings(phylotate::read_annotated(
					paste0(fasta,"_igphyml_jointpars_",run_id,".nex"),
					format=format))
			}
		}else{
			tf = suppressWarnings(read_annotated(paste0(fasta,append),format=format))
		}
		if(!is.null(states)){
			tf = condenseTrees(tf,states,palette)
		}
		germ = tf$tip.label[grep("_GERM",tf$tip.label)]
		tf$name = strsplit(germ,split="_")[[1]][1]
		tf = rerootGermline(tf,germ)
		trees[[i-1]] = tf
	}
	return(trees)
}

#write bootstapped lineage file input into igphyml
writeBootstrapFile = function(data,trees,dir=".",id="",rep=NULL,trait=NULL){
	file = paste0(dir,"/",id,"_blineages_pars.tsv")
	if(!is.null(rep)){
		file = paste0(dir,"/",id,"_blineages_",rep,"_pars.tsv")
	}
	outdir = paste0(dir,"/",id,"_",rep)
	#print(paste("outdir: ",outdir))
	dir.create(dir,showWarnings=FALSE)
	dir.create(outdir,showWarnings=FALSE)
	write(length(data),file=file)
	for(i in 1:length(data)){
		tree = ape::multi2di(trees[[i]])
		germid = paste0(i,"_GERM")
		tree$tip.label[which(tree$tip.label == "Germline")] = germid
		fastafile = paste0(outdir,"/",i,".fasta")
		treefile = paste0(outdir,"/",i,".tree")
		f = writeFasta(data[[i]],fastafile,germid,trait,dummy=TRUE)
		ape::write.tree(tree,file=treefile)
		write(paste(fastafile,treefile,germid,"N",sep="\t"), file=file,
			append=TRUE)
	}
	return(file)
}


# Write lineage file for IgPhyML use
writeLineageFile = function(data,trees,dir=".",id="N",rep=NULL,trait=NULL){
	file = paste0(dir,"/",id,"_lineages_pars.tsv")
	if(!is.null(rep)){
		file = paste0(dir,"/",id,"_lineages_",rep,"_pars.tsv")
	}
	outdir = paste0(dir,"/",id,"_recon_",rep)
	dir.create(dir,showWarnings=FALSE)
	dir.create(outdir,showWarnings=FALSE)

	dnames = unlist(lapply(data,function(x)x@clone))
	tnames = unlist(lapply(trees,function(x)x$name))
	if(sum(tnames != dnames) != 0){
		trees = trees[order(match(tnames,dnames))]
	}
	write(length(data),file=file)
	for(i in 1:length(data)){
		tree = trees[[i]]
		fastafile = paste0(outdir,"/",data[[i]]@clone,".fasta")
		treefile = paste0(outdir,"/",data[[i]]@clone,".tree")
		germid = paste0(data[[i]]@clone,"_GERM")
		f = writeFasta(data[[i]],fastafile,germid,trait,dummy=TRUE)
		tree$tip.label[which(tree$tip.label == "Germline")] = germid
		tree = ape::multi2di(tree)
		ape::write.tree(tree,file=treefile)
		write(paste(fastafile,treefile,germid,"N",sep="\t"), file=file,
			append=TRUE)
	}
	return(file)
}

#wrapper for build phylip lineage
buildPhylo = function(clone,trait,dnapars,temp_path=NULL,verbose=FALSE,rm_temp=TRUE,seq="SEQUENCE"){
	#stop("Dowser is not yet compatible with buildPhylipLineage.")
	if(seq != "SEQUENCE"){
		clone@data$SEQUENCE = clone@data[[seq]]
		if(seq == "HLSEQUENCE"){
		clone@germline = clone@hlgermline
		}else if(seq == "LSEQUENCE"){
			clone@germline = clone@lgermline
		}
	}
	phylo = tryCatch({
		alakazam::buildPhylipLineage(clone,dnapars,rm_temp=rm_temp,phylo=TRUE,verbose=verbose,
			temp_path=temp_path)},
		error=function(e){print(paste("buildPhylipLineage error:",e));stop()})
	phylo$name = clone@clone
	phylo$tree_method = "phylip::dnapars"
	phylo$edge_type = "genetic_distance"
	phylo$seq = seq
	phylo
}

#wrapper for pratchet + acctran
buildPratchet = function(clone,seq="SEQUENCE"){
	seqs = clone@data[[seq]]
	names = clone@data$SEQUENCE_ID
	if(seq == "HLSEQUENCE"){
		germline = clone@hlgermline
	}else if(seq == "LSEQUENCE"){
		germline = clone@lgermline
	}else{
		germline = clone@germline
	}
	seqs = base::append(seqs,germline)
	names = c(names,"Germline")
	seqs = strsplit(seqs,split="")
	names(seqs) = names
	data = phangorn::phyDat(ape::as.DNAbin(t(as.matrix(dplyr::bind_rows(seqs)))))
	tree = tryCatch(phangorn::pratchet(data,trace=FALSE),warning=function(w)w)
	tree = phangorn::acctran(ape::multi2di(tree),data)
	tree = rerootGermline(tree,"Germline",resolve=TRUE)
	tree$edge.length = tree$edge.length/nchar(germline)
	tree$name = clone@clone
	tree$tree_method = "phangorn::prachet"
	tree$edge_type = "genetic_distance"
	tree$seq = seq
	return(tree)
}

#### Preprocessing functions ####

#' Generate an ordered list of ChangeoClone objects for lineage construction
#' 
#' \code{formatClones} takes a \code{data.frame} or \code{tibble} with AIRR or Change-O style columns as input and 
#' masks gap positions, masks ragged ends, removes duplicates sequences, and merges 
#' annotations associated with duplicate sequences. If specified, it will un-merge duplicate sequences with different values specified in the \code{trait} option.
#' It returns a list of \code{ChangeoClone} objects ordered by number of sequences which serve as input for lineage reconstruction.
#' 
#' @param    data         data.frame containing the AIRR or Change-O data for a clone. See Details
#'                        for the list of required columns and their default values.
#' @param    id           name of the column containing sequence identifiers.
#' @param    seq          name of the column containing observed DNA sequences. All 
#'                        sequences in this column must be multiple aligned.
#' @param    germ         name of the column containing germline DNA sequences. All entries 
#'                        in this column should be identical for any given clone, and they
#'                        must be multiple aligned with the data in the \code{seq} column.
#' @param    vcall        name of the column containing V-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    jcall        name of the column containing J-segment allele assignments. All 
#'                        entries in this column should be identical to the gene level.
#' @param    junc_len     name of the column containing the length of the junction as a 
#'                        numeric value. All entries in this column should be identical 
#'                        for any given clone.
#' @param    clone        name of the column containing the identifier for the clone. All 
#'                        entries in this column should be identical.
#' @param    mask_char    character to use for masking and padding.
#' @param    max_mask     maximum number of characters to mask at the leading and trailing
#'                        sequence ends. If \code{NULL} then the upper masking bound will 
#'                        be automatically determined from the maximum number of observed 
#'                        leading or trailing Ns amongst all sequences. If set to \code{0} 
#'                        (default) then masking will not be performed.
#' @param    pad_end      if \code{TRUE} pad the end of each sequence with \code{mask_char}
#'                        to make every sequence the same length.
#' @param    text_fields  text annotation columns to retain and merge during duplicate removal.
#' @param    num_fields   numeric annotation columns to retain and sum during duplicate removal.
#' @param    seq_fields   sequence annotation columns to retain and collapse during duplicate 
#'                        removal. Note, this is distinct from the \code{seq} and \code{germ} 
#'                        arguments, which contain the primary sequence data for the clone
#'                        and should not be repeated in this argument.
#' @param    add_count    if \code{TRUE} add an additional annotation column called 
#'                        \code{COLLAPSE_COUNT} during duplicate removal that indicates the 
#'                        number of sequences that were collapsed.
#' @param    verbose      passed on to \code{collapseDuplicates}. If \code{TRUE}, report the 
#'                        numbers of input, discarded and output sequences; otherwise, process
#'                        sequences silently.
#' @param     collapse    iollapse identical sequences?
#' @param     traits      column ids to keep distinct during sequence collapse 
#' @param     region      if HL, include light chain information if available.
#' @param     heavy       name of heavy chain locus (default = "IGH")
#' @param     cell        name of the column containing cell assignment information
#' @param     locus       name of the column containing locus information
#' @param     nproc		  Number of cores to parallelize formating over.                        
#'
#' @return   A list of \link{ChangeoClone} objects containing modified clones.
#'
#' @details
#' This function is largely a wrapper for alakazam::makeChangeoClone.
#' The input data.frame (\code{data}) must columns for each of the required column name 
#' arguments: \code{id}, \code{seq}, \code{germ}, \code{vcall}, \code{jcall}, 
#' \code{junc_len}, and \code{clone}.  The default values are as follows:
#' \itemize{
#'   \item  \code{id       = "sequence_id"}:           unique sequence identifier.
#'   \item  \code{seq      = "sequence_alignment"}:         IMGT-gapped sample sequence.
#'   \item  \code{germ     = "germline_alignment_d_mask"}:  IMGT-gapped germline sequence.
#'   \item  \code{vcall    = "v_call"}:                V-segment allele call.
#'   \item  \code{jcall    = "j_call"}:                J-segment allele call.
#'   \item  \code{junc_len = "junction_length"}:       junction sequence length.
#'   \item  \code{clone    = "clone_id"}:                 clone identifier.
#' }
#' Additional annotation columns specified in the \code{text_fields}, \code{num_fields} 
#' or \code{seq_fields} arguments will be retained in the \code{data} slot of the return 
#' object, but are not required. If the input data.frame \code{data} already contains a 
#' column named \code{SEQUENCE}, which is not used as the \code{seq} argument, then that 
#' column will not be retained.
#' 
#' The default columns are IMGT-gapped sequence columns, but this is not a requirement. 
#' However, all sequences (both observed and germline) must be multiple aligned using
#' some scheme for both proper duplicate removal and lineage reconstruction. 
#'
#' The value for the germline sequence, V-segment gene call, J-segment gene call, 
#' junction length, and clone identifier are determined from the first entry in the 
#' \code{germ}, \code{vcall}, \code{jcall}, \code{junc_len} and \code{clone} columns, 
#' respectively. For any given clone, each value in these columns should be identical.
#'  
#' @seealso  Executes in order \code{alakazam::maskSeqGaps}, \code{alakazam::maskSeqEnds}, 
#'           \code{alakazam::padSeqEnds}, \code{alakazam::collapseDuplicates},
#' 			 and \code{processClones}. Returns a list of \link{ChangeoClone} objects 
#' 			which serve as input to \link{getTrees} and \link{bootstrapTrees}.
#' 
#' @examples
#' \dontrun{
#' data(ExampleDb)
#' clones = formatClones(ExampleDb,trait="sample_id")
#' }
#' @export
formatClones <- function(data, id="sequence_id", seq="sequence_alignment", 
                             germ="germline_alignment_d_mask", vcall="v_call", jcall="j_call",
                             junc_len="junction_length", clone="clone_id", mask_char="N",
                             max_mask=0, pad_end=FALSE, text_fields=NULL, num_fields=NULL, seq_fields=NULL,
                             add_count=TRUE, verbose=FALSE, nproc=1, collapse=TRUE,
                             region="H", heavy=NULL, cell="cell", locus="locus", traits=NULL) {

    clones <- data %>%
        dplyr::group_by(!!rlang::sym(clone)) %>%
        dplyr::do(DATA=alakazam::makeChangeoClone(.,
            id=id, seq=seq, germ=germ, vcall=vcall, jcall=jcall, junc_len=junc_len,
            clone=clone, mask_char=mask_char, max_mask=max_mask, pad_end=pad_end,
            text_fields=text_fields, num_fields=num_fields, seq_fields=seq_fields,
            add_count=add_count, verbose=verbose,collapse=collapse,
            region=region, heavy=heavy, cell=cell, locus=locus, traits=traits))

    if(region == "HL"){
    	seq_name = "HLSEQUENCE"
    }else{
    	seq_name = "SEQUENCE"
    }
    
    fclones = processClones(clones, nproc=nproc, seq=seq_name)
    fclones
}


# Clean clonal alignments and deduplicate based on trait
processClones = function(clones,nproc=1,seq){
	if(!"tbl" %in% class(clones)){
		print(paste("clones is of class",class(clones)))
		stop("clones must be a tibble of ChangeoClone objects!")
	}else{
		if(class(clones$DATA[[1]]) != "ChangeoClone"){
			print(paste("clones is list of class",class(clones$DATA[[1]])))
			stop("clones$DATA must be a list of ChangeoClone objects!")
		}
	}

	threshold = unlist(lapply(clones$DATA,function(x)
	length(x@data[[seq]]) >= 2))
	clones = clones[threshold,]

	clones$DATA = lapply(clones$DATA,function(x){
		x@data$SEQUENCE_ID=	gsub(":","_",x@data$SEQUENCE_ID);
		x
		})
	clones$DATA = lapply(clones$DATA,function(x){
		x@data$SEQUENCE_ID=gsub(";","_",x@data$SEQUENCE_ID);
		x
		})
	
	clones$DATA = lapply(clones$DATA,function(x){
		x@data$SEQUENCE_ID=gsub(",","_",x@data$SEQUENCE_ID);
		x
		})

	max = max(unlist(lapply(clones$DATA,function(x)max(nchar(x@data$SEQUENCE_ID)))))
	if(max > 1000){
		wc = which.max(unlist(lapply(clones$DATA,function(x)
			max(nchar(x@data$SEQUENCE_ID)))))
		stop(paste("Sequence ID of clone",clones$DATA[[wc]]@clone,"index",
			wc,"too long - over 1000 characters!"))
	}

	or = order(unlist(lapply(clones$DATA,function(x) nrow(x@data))),
		decreasing=TRUE)
	clones = clones[or,]

	clones$DATA = parallel::mclapply(clones$DATA,
		function(x)cleanAlignment(x,seq),mc.cores=nproc)

	if(.hasSlot(clones$DATA[[1]],"locus")){
		clones$LOCUS = unlist(lapply(clones$DATA,function(x)paste(x@locus,collapse=",")))
	}
	clones$SEQS = unlist(lapply(clones$DATA,function(x)nrow(x@data)))
	clones = dplyr::rowwise(clones)
	clones = dplyr::ungroup(clones)
	clones
}

# Reroot phylogenetic tree to have its germline sequence at a zero-length branch 
# to a node which is the direct ancestor of the tree's UCA. Assigns \code{uca}
# to be the ancestral node to the tree's germline sequence, as \code{germid} as
# the tree's germline sequence ID. Lifted from alakazam.
#
# @param   tree     An ape \code{phylo} object
# @param   germid   ID of the tree's predicted germline sequence
# @param   resolve  If \code{TRUE} reroots tree to specified germline sequnece.
#                   usually not necessary with IgPhyML trees analyzed with HLP model.
rerootGermline <- function(tree, germid, resolve=FALSE){
    if(resolve) {
        tree <- ape::root(tree, outgroup=germid, resolve.root=T, edge.label=TRUE)
    }
    tree <- ape::reorder.phylo(tree, "postorder")  
    edges <- tree$edge
    rootnode <- which(tree$tip.label==germid)
    rootedge <- which(edges[, 2] == rootnode)
    rootanc <- edges[edges[, 2] == rootnode, 1]
    mrcaedge <- which(edges[, 1] == rootanc & edges[, 2] != rootnode)
    if(length(mrcaedge) > 1){
            print("POLYTOMY AT ROOT?!")
            quit(save="no", status=1, runLast=FALSE)
    }
    tree$edge.length[mrcaedge] <- tree$edge.length[mrcaedge] + tree$edge.length[rootedge]
    tree$edge.length[rootedge] <- 0
    tree$uca <- rootanc
    tree$germid <- germid
    
    return(tree)
}

#' Estimate lineage tree topologies, branch lengths,
#' and internal node states if desired
#' 
#' \code{getTrees} Tree building function.
#' @param    clones  	a tiblle of \code{changeoClone} objects, the output of \link{formatClones}
#' @param    data  		list of \code{changeoClone} objects, alternate input
#' @param    trait 		trait to use for parsimony models (required if \code{igphyml} specified)
#' @param 	 build	    program to use for tree building (phangorn, dnapars)
#' @param 	 exec	    location of desired phylogenetic executable
#' @param    igphyml 	location of igphyml executible if trait models desired (optional)
#' @param    id 		unique identifer for this analysis (required if \code{igphyml} or \code{dnapars} specified)
#' @param    dir    	directory where temporary files will be placed (required if \code{igphyml} or \code{dnapars} specified)
#' @param    modelfile 	file specifying parsimony model to use
#' @param    trees 		tree topologies to use if laready available
#' @param    nproc 		number of cores to parallelize computations
#' @param    quiet		amount of rubbish to print to console
#' @param    rm_temp	remove temporary files (default=TRUE)
#' @param    palette 	a named vector specifying colors for each state
#' @param    resolve 	how should polytomies be resolved?
#' @param    seq        column name containing sequence information
#' 
#' @return   A list of \code{phylo} objects in the same order as \code{data}.
#'
#' @details
#' Estimates phylogenetic tree topologies and branch lengths for a list of \code{changeoClone} objects.
#' By default, it will use phangnorn::pratchet to estimate maximum parsimony tree topologies, and 
#' ape::acctran to estimate branch lengths. If \code{igpyhml} is specified, internal node \code{trait} values
#' will be predicted by maximum parsimony. In this case, \code{dir} will need to be specified as a temporary
#' directory to place all the intermediate files (will be created if not available). Further, \code{id} will
#' need to specified to serve as a unique identifier for the temporary files. This should be chosen to ensure
#' that multiple \code{getTrees} calls using the same \code{dir} do not overwrite each others files. 
#' 
#' \code{modelfile} is written automatically if not specified, but doesn't include any constraints. Intermediate
#' files are deleted by default. This can be toggled using (\code{rm_files}).
#'  
#' @seealso \link{formatClones}, \link{bootstrapTrees}
#' @examples
#' \dontrun{
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb,trait="sample_id")
#'
#' trees = getTrees(clones[1:2])
#' plotTrees(trees[[1]])
#' 
#' trees = getTrees(clones[1:2],igphyml=igphyml,id="temp",dir="temp",trait="sample_id")
#' plotTrees(trees[[1]])
#' }
#' @export
getTrees = function(clones,data=NULL,trait=NULL,id=NULL,dir=NULL,modelfile=NULL,
	build="pratchet",exec=NULL,igphyml=NULL,trees=NULL,nproc=1,quiet=0,rm_temp=TRUE,
	palette=NULL,resolve=2,seq=NULL){

	data = clones$DATA
	if(is.null(id)){
            id <- "sample"
    }
    if(class(data) != "list"){
        data <- list(data)
    }
    if(class(data[[1]]) != "ChangeoClone"){
        stop("Input data must be a list of ChangeoClone objects")
    }
    big <- FALSE
    if(sum(unlist(lapply(data, function(x)nrow(x@data)))) > 10000){
        big <- TRUE
    }
    if(!rm_temp && big){
        warning("Large dataset - best to set rm_temp=TRUE")
    }
	if(!is.null(igphyml)){
		igphyml <- path.expand(igphyml)
        if(!is.null(dir)){
            if(!dir.exists(dir)){
                dir.create(dir)
            }
        }else{
            dir <- alakazam::makeTempDir(id)
            if(big){
                warning("Large dataset - best to set dir and id params")
            }
        }
        if(is.null(trait)){
            stop("trait must be specified when running igphyml")
        }
		if(is.null(modelfile)){
			states = unique(unlist(lapply(data,function(x)x@data[,trait])))
			modelfile = makeModelFile(states,file=paste0(dir,"/",id,"_modelfile.txt"))
		}else{
			states = readModelFile(modelfile)
		}
		#if igphyml is specified, append trait value to sequence ids
		data = lapply(data,function(x){
			x@data$SEQUENCE_ID = paste0(x@data$SEQUENCE_ID,"_",x@data[[trait]])
			x})
	}
	if(build=="dnapars"){
		if(!is.null(dir)){
            if(!dir.exists(dir)){
                dir.create(dir)
            }
        }else{
            dir <- alakazam::makeTempDir(id)
            if(big){
                warning("Large dataset - best to set dir and id params")
            }
        }
	}

	if(class(data) != "list"){
		data = list(data)
	}
	if(!is.null(dir)){
		if(!dir.exists(dir)){
			dir.create(dir)
		}
	}
	rm_dir = NULL
	if(rm_temp){
		rm_dir=paste0(dir,"/",id,"_recon_trees")
	}
	
	if(is.null(trees)){
		reps = as.list(1:length(data))
		if(is.null(seq)){
			seqs = unlist(lapply(data,function(x)x@phylo_seq))
		}else{
			seqs = rep(seq,length=length(data))
		}
		if(build=="dnapars"){
			trees = parallel::mclapply(reps,function(x)
				buildPhylo(data[[x]],
					trait,exec,
					temp_path=paste0(dir,"/",id,"_trees_",x),
					rm_temp=rm_temp,
					seq=seqs[x]),
				mc.cores=nproc)
		}else{
			trees = parallel::mclapply(reps,function(x)
				buildPratchet(data[[x]],seq=seqs[x]),
				mc.cores=nproc)
		}
	}
	
	if(!is.null(igphyml)){
		file = writeLineageFile(data=data, trees=trees, dir=dir,
			id=id, trait=trait, rep="trees")
	
		mtrees = reconIgPhyML(file, modelfile, igphyml=igphyml, 
			mode="trees", cloneid=NULL, quiet=quiet, nproc=nproc,
			rm_files=rm_temp, rm_dir=rm_dir, states=states, 
			palette=palette,resolve=resolve)
	}else{
		mtrees = trees
	}
	clones$TREE = mtrees
	clones
}



#' Estimate lineage tree topologies, branch lengths,
#' and internal node states if desired
#' 
#' \code{getTrees} Tree building function.
#' @param    clones  	a tiblle of \code{changeoClone} objects, the output of \link{formatClones}
#' @param    edge_type  Either \code{genetic_distance} (mutations per site) or \code{mutations}   
#' 
#' @return   A tibble with \code{phylo} objects that have had branch lengths rescaled.
#'
#' @details
#' Uses clones$TREE[[1]]$edge_type to determine how branches are currently scaled.
#'  
#' @seealso \link{getTrees}
#' @export
scaleBranches = function(clones,edge_type="mutations"){
	if(!"tbl" %in% class(clones)){
		print(paste("clones is of class",class(clones)))
		stop("clones must be a tibble of ChangeoClone objects!")
	}else{
		if(class(clones$DATA[[1]]) != "ChangeoClone"){
			print(paste("clones is list of class",class(clones$DATA[[1]])))
			stop("clones must be a list of ChangeoClone objects!")
		}
	}
	if(!"TREE" %in% names(clones)){
		stop("clones must have TREE column!")
	}
	# Need to add a variable storing whether branch lengths equal mutations or genetic distance
	lengths = unlist(lapply(1:length(clones$TREE),
		function(x){
		if(clones$TREE[[x]]$seq == "HLSEQUENCE"){
			return(nchar(clones$DATA[[x]]@hlgermline))
		}else if(clones$TREE[[x]]$seq == "LSEQUENCE"){
			return(nchar(clones$DATA[[x]]@lgermline))
		}else{
			return(nchar(clones$DATA[[x]]@germline))
		}}))

	TREE = lapply(1:length(clones$TREE),function(x){
		if(clones$TREE[[x]]$edge_type == "mutations" && edge_type == "genetic_distance"){
			clones$TREE[[x]]$edge.length = clones$TREE[[x]]$edge.length/lengths[x]
			clones$TREE[[x]]$edge_type = "genetic_distance"
			clones$TREE[[x]]
		}else if(clones$TREE[[x]]$edge_type == "genetic_distance" && edge_type == "mutations"){
			clones$TREE[[x]]$edge.length = clones$TREE[[x]]$edge.length*lengths[x]
			clones$TREE[[x]]$edge_type = "mutations"
			clones$TREE[[x]]
		}else{
			clones$TREE[[x]]
		}})
			
	clones$TREE = TREE
	clones
}



#' Create a bootstrap distribution for clone sequence alignments, and estimate trees for eac
#' bootstrap replicate.
#' 
#' \code{bootstrapTrees} Phylogenetic bootstrap function.
#' @param    clones  	tibble \code{changeoClone} objects, the output of \link{formatClones}
#' @param    bootstraps number of bootstrap replicates to perform
#' @param    trait 		trait to use for parsimony models (required if \code{igphyml} specified)
#' @param 	 build	    program to use for tree building (phangorn, dnapars)
#' @param 	 exec	    location of desired phylogenetic executable
#' @param    igphyml 	location of igphyml executible if trait models desired (optional)
#' @param    id 		unique identifer for this analysis (required if \code{igphyml} or \code{dnapars} specified)
#' @param    dir    	directory where temporary files will be placed (required if \code{igphyml} or \code{dnapars} specified)
#' @param    modelfile 	file specifying parsimony model to use
#' @param    trees 		tree topologies to use if aready available (bootstrapping will not be perfomed)
#' @param    nproc 		number of cores to parallelize computations
#' @param    quiet		amount of rubbish to print to console
#' @param    rm_temp	remove temporary files (default=TRUE)
#' @param    palette 	a named vector specifying colors for each state
#' @param    resolve 	how should polytomies be resolved?
#' @param    keeptrees  keep trees estimated from bootstrap replicates? (TRUE)
#' @param    lfile      lineage file input to igphyml if desired (experimental)
#' @param    rep  		current bootstrap replicate (experimental)
#' @param    seq        column name containing sequence information
#'
#' @return   A list of trees and/or switch counts for each bootstrap replicate.
#'
#' @details
#' Tree building details are the same as \link{getTrees}. 
#' If \code{keeptrees=TRUE} (default) the returned object will contain a list named "trees"
#' which contains a list of estimated tree objects for each bootstrap replicate. The object is
#' structured like: trees[[<replicate>]][[<tree index>]].
#' If \code{igphyml} is specified (as well as \code{trait}, \code{temp}, and \code{id}), the
#' returned object will contain a \code{tibble} named "switches" containing switch count information.
#' This object can be passed to \link{PStest} and other functions to perform parsimony based trait value
#' tests. 
#'  
#' @seealso Uses output from \link{formatClones} with similar arguments to \link{getTrees}. Output can be 
#' visualized with \link{plotTrees}, and tested with \link{PStest}, \link{SCtest}, and \link{SPtest}.
#' 
#' @examples
#' \dontrun{
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb,trait="sample_id")
#' 
#' btrees = bootstrapTrees(clones[1:2],bootstraps=100)
#' plotTrees(btrees$trees[[4]][[1]])
#' 
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' btrees = bootstrapTrees(clones[1:2],bootstraps=100,nproc=1,
#'	igphyml=igphyml,trait="sample_id",id="temp",dir="temp")
#' plotTrees(btrees$trees[[4]][[1]])
#' PStest(btrees$switches)
#' }
#' @export
bootstrapTrees = function(clones, bootstraps, nproc=1, trait=NULL, dir=NULL, 
	id=NULL, modelfile=NULL,build="pratchet",exec=NULL,igphyml=NULL,trees=NULL,
	quiet=0,rm_temp=TRUE,palette=NULL,resolve=2,rep=NULL,
	keeptrees=TRUE, lfile=NULL, seq="SEQUENCE"){
	data = clones$DATA
	if(is.null(id)){
            id <- "sample"
    }
	if(class(data) != "list"){
		data = list(data)
	}
	if(!is.null(dir)){
		if(!dir.exists(dir)){
			dir.create(dir)
		}
	}
	if(class(data[[1]]) != "ChangeoClone"){
        stop("Input data must be a list of ChangeoClone objects")
    }
    big <- FALSE
    if(sum(unlist(lapply(data, function(x)nrow(x@data)))) > 10000){
        big <- TRUE
    }
    if(!rm_temp && big){
        warning("Large dataset - best to set rm_temp=TRUE")
    }
	if(!is.null(igphyml)){
		igphyml = path.expand(igphyml)
		if(!is.null(dir)){
            if(!dir.exists(dir)){
                dir.create(dir)
            }
        }else{
            dir <- alakazam::makeTempDir(id)
            if(big){
                warning("Large dataset - best to set dir and id params")
            }
        }
        if(is.null(trait)){
            stop("trait must be specified when running igphyml")
        }
		if(is.null(modelfile)){
			states = unique(unlist(lapply(data,function(x)x@data[,trait])))
			modelfile = makeModelFile(states,file=paste0(dir,"/",id,"_modelfile.txt"))
		}else{
			states = readModelFile(modelfile)
		}
		#if igphyml is specified, append trait value to sequence ids
		data = lapply(data,function(x){
			x@data$SEQUENCE_ID = paste0(x@data$SEQUENCE_ID,"_",x@data[[trait]])
			x})
	}
	if(build=="dnapars"){
		if(is.null(dir) || is.null(id)){
			stop("dir, and id parameters must be specified when running dnapars")
		}
	}
	if(is.null(rep)){
		reps = as.list(1:bootstraps)
		l = parallel::mclapply(reps,function(x)
			bootstrapTrees(clones,rep=x, 
			trait=trait, modelfile=modelfile,build=build, 
			exec=exec, igphyml=igphyml, 
			id=id, dir=dir, bootstraps=bootstraps,
			nproc=1, rm_temp=rm_temp, quiet=quiet,
			trees=trees,resolve=resolve,keeptrees=keeptrees,
			lfile=lfile,seq=seq),
			mc.cores=nproc)
		results = list()
		results$switches = NULL
		results$trees = NULL
		if(!is.null(igphyml)){
			results$switches = dplyr::bind_rows(lapply(l,function(x)x$switches))
		}
		if(keeptrees){
			results$trees = lapply(l,function(x)x$trees)
		}
		if(rm_temp){
			if(file.exists(paste0(dir,"/",id,"_modelfile.txt"))){
				unlink(paste0(dir,"/",id,"_modelfile.txt"))
			}
		}
		return(results)
	}else{
		rm_dir=paste0(dir,"/",id,"_recon_",rep)
		if(is.null(trees)){
			for(i in 1:length(data)){
				if(quiet > 3){
					print(table(data[[i]]@data[,trait]))
				}	
				data[[i]] = bootstrapClones(data[[i]], reps=1)[[1]]
			}
			if(quiet > 1){print("building trees")}
			reps = as.list(1:length(data))
			if(is.null(seq)){
			      seqs = unlist(lapply(data,function(x)x@phylo_seq))
		    }else{
		          seqs = rep(seq,length=length(data))
		    }
			if(build=="dnapars"){
				trees = lapply(reps,function(x)
					buildPhylo(data[[x]],
						trait,exec,
						temp_path=paste0(dir,"/",id,"_",rep,"_trees_",x),
						rm_temp=rm_temp,seq=seq))
			}else{
				trees = parallel::mclapply(reps,function(x)
					buildPratchet(data[[x]],seq),
					mc.cores=nproc)
			}
		}
		results = list()
		if(!is.null(igphyml)){
			if(is.null(lfile)){	
				file = writeLineageFile(data=data, trees=trees, dir=dir,
					id=id, trait=trait, rep=rep)
			}else{
				file = lfile
			}
			rseed = floor(stats::runif(1,0,10^7))+rep
			switches = reconIgPhyML(file, modelfile, igphyml=igphyml, 
				mode="switches", type="recon", cloneid=rep, 
				quiet=quiet, rm_files=FALSE, rm_dir=NULL, nproc=nproc,
				resolve=resolve, rseed=rseed)
			permuted = reconIgPhyML(file, modelfile, igphyml=igphyml, 
				mode="switches", type="permute", cloneid=rep, 
				quiet=quiet, rm_files=FALSE, rm_dir=NULL, nproc=nproc,
				resolve=resolve, rseed=rseed)
			if(!rm_temp){rm_dir=NULL}
			permuteAll = reconIgPhyML(file, modelfile, igphyml=igphyml, 
				mode="switches", type="permuteAll", cloneid=rep, 
				quiet=quiet, rm_files=rm_temp, rm_dir=rm_dir, nproc=nproc,
				resolve=resolve, rseed=rseed)
			switches = rbind(switches,permuted)
			switches = rbind(switches,permuteAll)
			switches$ID = id
			temp = switches$REP
			switches$REP = switches$CLONE
			switches$CLONE = temp
			clone_index = unlist(lapply(data,function(x)x@clone))
			switches$CLONE = clone_index[switches$CLONE+1]
			results$switches = switches
		}
		if(keeptrees){
			results$trees = trees
		}
		return(results)
	}
}

