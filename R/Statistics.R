# Summary statistics and data transformation for switch count distributions

#' Performs PS (parsimony score) test on switch data
#' 
#' \code{PStest} performs a PS test
#' @param    switches     Data frame from bootstrapTrees
#' @param    bylineage    Perform test for each lineage individually? (FALSE)
#' @param    pseudocount  Pseudocount for P value calculations
#'
#' @return   A list containing a \code{tibble} with mean PS statistics, and another 
#' with PS statistics per repetition.
#' @details
#' Output data table columns:
#' RECON = PS for observed data
#' PERMUTE = PS for permuted data
#' DELTA = RECON - PERMUTE
#' PLT = p value for DELTA < 0
#' PGT = p value for DELTA < 0
#' \itemize{
#'   \item  \code{RECON}: PS for observed data.
#'   \item  \code{PERMUTE}: PS for permuted data.
#'   \item  \code{DELTA}:  RECON - PERMUTE.
#'   \item  \code{PLT}: p value that DELTA < 0
#'   \item  \code{PGT}: p value that DELTA > 0
#'   \item  \code{STAT}: Statistic used (PS).
#'   \item  \code{REP}: Bootstrap repetition.
#'   \item  \code{REPS}: Total number of ootstrap repetition.
#'}
#'  
#' @seealso Uses output from \link{bootstrapTrees}. Related to \link{SPtest}
#' and \link{SCtest}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb,trait="sample_id")
#' btrees = bootstrapTrees(clones[1:2],bootstraps=100,nproc=1,
#'	igphyml=igphyml,trait="sample_id",id="temp",dir="temp")
#' PStest(btrees$switches)
#' }
#' @export
PStest <- function(switches,bylineage=FALSE,pseudocount=0){
	switches <- switches %>% 
		dplyr::filter(!!rlang::sym("TO") != "N" & 
			!!rlang::sym("TO") != !!rlang::sym("FROM") &
			 !!rlang::sym("FROM") != "UCA")
	if(!bylineage){
		reps <- switches  %>%
			dplyr::group_by(!!rlang::sym("REP"), !!rlang::sym("TYPE")) %>% 
			dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES"))) %>% 
			dplyr::ungroup() %>% 
			tidyr::spread(!!rlang::sym("TYPE"),!!rlang::sym("SWITCHES")) %>%
			dplyr::mutate(DELTA = !!rlang::sym("RECON") - !!rlang::sym("PERMUTE"))
	}else{
		reps <- switches  %>%
			dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("CLONE"), !!rlang::sym("REP")) %>% 
			dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES"))) %>% 
			dplyr::ungroup() %>% 
			tidyr::spread(!!rlang::sym("TYPE"),!!rlang::sym("SWITCHES")) %>%
			dplyr::mutate(DELTA = !!rlang::sym("RECON") - !!rlang::sym("PERMUTE"))
	}

	if(!bylineage){
		means <- reps %>% 
			dplyr::summarize(
				RECON = mean(!!rlang::sym("RECON")),
				PERMUTE = mean(!!rlang::sym("PERMUTE")),
				PLT = (sum(!!rlang::sym("DELTA") >= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				PGT = (sum(!!rlang::sym("DELTA") <= 0) + pseudocount)/
					(dplyr::n() + pseudocount),	
				DELTA = mean(!!rlang::sym("DELTA")))	
	}else{
		means <- reps %>% 
			dplyr::group_by(!!rlang::sym("CLONE")) %>% 
			dplyr::summarize(
				RECON = mean(!!rlang::sym("RECON")),
				PERMUTE = mean(!!rlang::sym("PERMUTE")),
				PLT = (sum(!!rlang::sym("DELTA") >= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				PGT = (sum(!!rlang::sym("DELTA") <= 0) + pseudocount)/
					(dplyr::n() + pseudocount),	
				DELTA = mean(!!rlang::sym("DELTA")))	
	}

	means$STAT <- "PS"
	reps$STAT <- "PS"
	means$REPS = length(unique(reps$REP))
	results <- list()
	results$reps <- reps
	results$means <- means
	results
}

#' Performs SP (switch proportion) test on switch data
#' 
#' \code{SPtest} performs an SP test
#' @param    switches     Data frame from bootstrapTrees
#' @param    permuteAll   Permute among trees?
#' @param    from         Include only switches from this state?
#' @param    to           Include only switches to this state?
#' @param    dropzeros    Drop switches with zero counts?
#' @param    bylineage    Perform test for each lineage individually?
#' @param    pseudocount  Pseudocount for P value calculations
#'
#' @return   A list containing a \code{tibble} with mean SP statistics, and another 
#' with SP statistics per repetition.
#'
#' @details
#' Output data table columns:
#' RECON = SP for observed data
#' PERMUTE = SP for permuted data
#' DELTA = RECON - PERMUTE
#' PLT = p value for DELTA < 0
#' PGT = p value for DELTA < 0
#' \itemize{
#'   \item  \code{FROM}: State going from.
#'   \item  \code{TO}: State going to.
#'   \item  \code{RECON}: SP for observed data.
#'   \item  \code{PERMUTE}: SP for permuted data.
#'   \item  \code{DELTA}:  RECON - PERMUTE.
#'   \item  \code{PLT}: p value that DELTA < 0
#'   \item  \code{PGT}: p value that DELTA > 0
#'   \item  \code{STAT}: Statistic used (SP).
#'   \item  \code{REP}: Bootstrap repetition.
#'   \item  \code{REPS}: Total number of ootstrap repetition.
#'}
#'  
#' @seealso Uses output from \link{bootstrapTrees}. Related to \link{PStest}
#' and \link{SCtest}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb,trait="sample_id")
#' btrees = bootstrapTrees(clones[1:2],bootstraps=100,nproc=1,
#'	igphyml=igphyml,trait="sample_id",id="temp",dir="temp")
#' SPtest(btrees$switches)
#' }
#' @export
SPtest <- function(switches, permuteAll=FALSE, 
	from=NULL, to=NULL, dropzeros=TRUE,
	bylineage=FALSE, pseudocount=0){

	permute <- dplyr::quo(!!rlang::sym("PERMUTE"))
	if(permuteAll){
		permute <- dplyr::quo(!!rlang::sym("PERMUTEALL"))
	}

	switches <- switches %>% 
		dplyr::filter(!!rlang::sym("TO") != !!rlang::sym("FROM") & !!rlang::sym("FROM") != "UCA")

	if(!is.null(from)){
		from <- dplyr::enquo(from)
		switches <- dplyr::filter(switches,!!rlang::sym("FROM") == !!from)
	}
	if(!is.null(to)){
		to <- dplyr::enquo(to)
		switches <- dplyr::filter(switches,!!rlang::sym("TO") == !!to)
	}

	if(!bylineage){
		reps <- switches  %>%
			dplyr::group_by(!!rlang::sym("REP"), !!rlang::sym("TYPE"), !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
			dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES")))
		reps <- reps %>% 
			dplyr::filter(!!rlang::sym("TO") != "N") %>%
			dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("REP")) %>% 
			dplyr::mutate(PROP = !!rlang::sym("SWITCHES")/sum(!!rlang::sym("SWITCHES"))) %>%
			dplyr::select(-!!rlang::sym("SWITCHES")) %>% tidyr::spread(!!rlang::sym("TYPE"),!!rlang::sym("PROP"))
	}else{		
		reps <- switches  %>%
			dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("CLONE"), !!rlang::sym("REP")) %>% 
			dplyr::mutate(PROP = !!rlang::sym("SWITCHES")/sum(!!rlang::sym("SWITCHES"))) %>% 
			dplyr::select(-!!rlang::sym("SWITCHES")) %>% tidyr::spread(!!rlang::sym("TYPE"),!!rlang::sym("PROP"))
	}
	
	reps <- reps %>% dplyr::mutate(DELTA = !!rlang::sym("RECON")-!!permute)
	
	if(dropzeros){
		from_type <- 
			reps %>% dplyr::group_by(!!rlang::sym("FROM")) %>%
			dplyr::summarize(d = sum(!!rlang::sym("RECON"))) %>%
			dplyr::filter(!!rlang::sym("d") == 0) %>%
			.$FROM
		
		to_type <- 
			reps %>% dplyr::group_by(!!rlang::sym("TO")) %>%
			dplyr::summarize(d = sum(!!rlang::sym("RECON"))) %>%
			dplyr::filter(!!rlang::sym("d") == 0) %>%
			.$TO
		
		remove <- to_type[to_type %in% from_type]
		reps <- reps[!reps$FROM %in% remove &
			!reps$TO %in% remove,]
	}

	if(!bylineage){
		means <- reps %>% 
			dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
			dplyr::summarize(
				RECON = mean(!!rlang::sym("RECON")),
				PERMUTE = mean(!!permute),
				PLT = (sum(!!rlang::sym("DELTA") >= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				PGT = (sum(!!rlang::sym("DELTA") <= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				DELTA = mean(!!rlang::sym("DELTA")))
	}else{
		means <- reps %>% 
			dplyr::group_by(!!rlang::sym("CLONE"), !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
			dplyr::summarize(
				RECON = mean(!!rlang::sym("RECON")),
				PERMUTE = mean(!!permute),
				PLT = (sum(!!rlang::sym("DELTA") >= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				PGT = (sum(!!rlang::sym("DELTA") <= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				DELTA = mean(!!rlang::sym("DELTA")))		
	}

	means$STAT <- "SP"
	reps$STAT <- "SP"
	means$REPS = length(unique(reps$REP))
	results <- list()
	results$reps <- reps
	results$means <- means
	results
}

#' Performs SC (switch count) test on switch data
#' 
#' \code{SCtest} performs an SC test
#' @param    switches     Data frame from bootstrapTrees
#' @param    permuteAll   Permute among trees?
#' @param    from         Include only switches from this state?
#' @param    to           Include only switches to this state?
#' @param    dropzeros    Drop switches with zero counts?
#' @param    bylineage    Perform test for each lineage individually?
#' @param    pseudocount  Pseudocount for P value calculations
#'
#' @return   A list containing a \code{tibble} with mean SC statistics, and another 
#' with SC statistics per repetition.
#'
#' @details
#' Output data table columns:
#' RECON = SC for observed data
#' PERMUTE = SC for permuted data
#' DELTA = RECON - PERMUTE
#' PLT = p value for DELTA < 0
#' PGT = p value for DELTA < 0
#' \itemize{
#'   \item  \code{FROM}: State going from.
#'   \item  \code{TO}: State going to.
#'   \item  \code{RECON}: SC for observed data.
#'   \item  \code{PERMUTE}: SC for permuted data.
#'   \item  \code{DELTA}:  RECON - PERMUTE.
#'   \item  \code{PLT}: p value that DELTA < 0
#'   \item  \code{PGT}: p value that DELTA > 0
#'   \item  \code{STAT}: Statistic used (SC).
#'   \item  \code{REP}: Bootstrap repetition.
#'   \item  \code{REPS}: Total number of ootstrap repetition.
#'}
#'  
#' @seealso Uses output from \link{bootstrapTrees}. Related to \link{PStest}
#' and \link{SPtest}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb,trait="sample_id")
#' btrees = bootstrapTrees(clones[1:2],bootstraps=100,nproc=1,
#'	igphyml=igphyml,trait="sample_id",id="temp",dir="temp")
#' SCtest(btrees$switches)
#' }
#' @export
SCtest <- function(switches,dropzeros=TRUE,
	bylineage=FALSE, pseudocount=0, from=NULL, to=NULL,
	permuteAll=FALSE){

	permute <-dplyr::quo(!!rlang::sym("PERMUTE"))
	if(permuteAll){
		permute <-dplyr::quo(!!rlang::sym("PERMUTEALL"))
	}

	switches <- switches %>% 
		dplyr::filter(!!rlang::sym("TO") != !!rlang::sym("FROM") & !!rlang::sym("FROM") != "UCA")

	if(!is.null(from)){
		from <-dplyr::enquo(from)
		switches <- dplyr::filter(switches,!!rlang::sym("FROM") == !!from)
	}
	if(!is.null(to)){
		to <-dplyr::enquo(to)
		switches <- dplyr::filter(switches,!!rlang::sym("TO") == !!to)
	}

	if(!bylineage){
		reps <- switches  %>%
			dplyr::group_by(!!rlang::sym("REP"), !!rlang::sym("TYPE"), !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
			dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES")))

		reps <- reps %>% 
			dplyr::filter(!!rlang::sym("TO") != "N") %>%
			dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("REP")) %>% 
			dplyr::mutate(COUNT = !!rlang::sym("SWITCHES")) %>%
			dplyr::select(-!!rlang::sym("SWITCHES")) %>% tidyr::spread(!!rlang::sym("TYPE"),!!rlang::sym("COUNT"))
	}else{
		reps <- switches  %>%
			dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("CLONE"), !!rlang::sym("REP")) %>% 
			dplyr::mutate(COUNT = !!rlang::sym("SWITCHES")) %>% 
			dplyr::select(-!!rlang::sym("SWITCHES")) %>% tidyr::spread(!!rlang::sym("TYPE"),!!rlang::sym("COUNT"))
	}

	reps <- reps %>% dplyr::mutate(DELTA = !!rlang::sym("RECON")-!!permute)

	if(dropzeros){
		from_type <- 
			reps %>% dplyr::group_by(!!rlang::sym("FROM")) %>%
			dplyr::summarize(d = sum(!!rlang::sym("RECON"))) %>%
			dplyr::filter(!!rlang::sym("d") == 0) %>%
			.$FROM
		
		to_type <- 
			reps %>% dplyr::group_by(!!rlang::sym("TO")) %>%
			dplyr::summarize(d = sum(!!rlang::sym("RECON"))) %>%
			dplyr::filter(!!rlang::sym("d") == 0) %>%
			.$TO
		
		remove <- to_type[to_type %in% from_type]
		reps <- reps[!reps$FROM %in% remove &
			!reps$TO %in% remove,]
	}

	if(!bylineage){
		means <- reps %>% 
			dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
			dplyr::summarize(
				RECON = mean(!!rlang::sym("RECON")),
				PERMUTE = mean(!!permute),
				PLT = (sum(!!rlang::sym("DELTA") >= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				PGT = (sum(!!rlang::sym("DELTA") <= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				DELTA = mean(!!rlang::sym("DELTA")))
	}else{
		means <- reps %>% 
			dplyr::group_by(!!rlang::sym("CLONE"), !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
			dplyr::summarize(
				RECON = mean(!!rlang::sym("RECON")),
				PERMUTE = mean(!!permute),
				PLT = (sum(!!rlang::sym("DELTA") >= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				PGT = (sum(!!rlang::sym("DELTA") <= 0) + pseudocount)/
					(dplyr::n() + pseudocount),
				DELTA = mean(!!rlang::sym("DELTA")))		
	}
	means$STAT <- "SC"
	reps$STAT <- "SC"
	means$REPS = length(unique(reps$REP))
	results <- list()
	results$reps <- reps
	results$means <- means
	results
}
