# Summary statistics and data transformation for switch count distributions

#' Performs PS (parsimony score) test on switch data
#' 
#' \code{testPS} performs a PS test
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
#' @seealso Uses output from \link{bootstrapTrees}. Related to \link{testSP}
#' and \link{testSC}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleDb)
#' ExampleDb$sample_id <- sample(ExampleDb$sample_id)
#' clones <- formatClones(ExampleDb, trait="sample_id")
#' btrees <- bootstrapTrees(clones[1:2], bootstraps=100, nproc=1,
#'    igphyml=igphyml, trait="sample_id", id="temp", dir="temp")
#' testPS(btrees$switches)
#' }
#' @export
testPS <- function(switches, bylineage=FALSE, pseudocount=0,
    alternative=c("less","two.sided","greater")){
    switches <- switches %>% 
        dplyr::filter(!!rlang::sym("TO") != "N" & 
            !!rlang::sym("TO") != !!rlang::sym("FROM") &
             !!rlang::sym("FROM") != "UCA")
    if(!bylineage){
        reps <- switches  %>%
            dplyr::group_by(!!rlang::sym("REP"), !!rlang::sym("TYPE")) %>% 
            dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES"))) %>% 
            dplyr::ungroup() %>% 
            tidyr::spread(!!rlang::sym("TYPE"), !!rlang::sym("SWITCHES")) %>%
            dplyr::mutate(DELTA = !!rlang::sym("RECON") - !!rlang::sym("PERMUTE"))
    }else{
        reps <- switches  %>%
            dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("CLONE"), 
            	!!rlang::sym("REP")) %>% 
            dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES"))) %>% 
            dplyr::ungroup() %>% 
            tidyr::spread(!!rlang::sym("TYPE"), !!rlang::sym("SWITCHES")) %>%
            dplyr::mutate(DELTA = !!rlang::sym("RECON") - !!rlang::sym("PERMUTE"))
    }

    if(!bylineage){
        if(alternative[1] == "two.sided"){
            means <- reps %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "greater"){
            means <- reps %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),
                    DELTA = mean(!!rlang::sym("DELTA")))
        }
    }else{
        if(alternative[1] == "two.sided"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("CLONE")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "greater"){
            means <- reps %>% 
                dplyr::group_by(!!rlang::sym("CLONE")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>% 
                dplyr::group_by(!!rlang::sym("CLONE")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),
                    DELTA = mean(!!rlang::sym("DELTA")))
        }            
    }

    means$STAT <- "PS"
    reps$STAT <- "PS"
    means$REPS <- length(unique(reps$REP))
    results <- list()
    results$reps <- reps
    results$means <- means
    results
}

#' Performs SP (switch proportion) test on switch data
#' 
#' \code{testSP} performs an SP test
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
#' @seealso Uses output from \link{bootstrapTrees}. Related to \link{testPS}
#' and \link{testSC}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb, trait="sample_id")
#' btrees = bootstrapTrees(clones[1:2], bootstraps=100, nproc=1,
#'    igphyml=igphyml, trait="sample_id", id="temp", dir="temp")
#' testSP(btrees$switches)
#' }
#' @export
testSP <- function(switches, permuteAll=FALSE, 
    from=NULL, to=NULL, dropzeros=TRUE,
    bylineage=FALSE, pseudocount=0, alternative=c("two.sided","greater","less")){

    permute <- dplyr::quo(!!rlang::sym("PERMUTE"))
    if(permuteAll){
        permute <- dplyr::quo(!!rlang::sym("PERMUTEALL"))
    }

    switches <- switches %>% 
        dplyr::filter(!!rlang::sym("TO") != !!rlang::sym("FROM") & 
        	!!rlang::sym("FROM") != "UCA")

    if(!is.null(from)){
        from <- dplyr::enquo(from)
        switches <- dplyr::filter(switches, !!rlang::sym("FROM") == !!from)
    }
    if(!is.null(to)){
        to <- dplyr::enquo(to)
        switches <- dplyr::filter(switches, !!rlang::sym("TO") == !!to)
    }

    if(!bylineage){
        reps <- switches  %>%
            dplyr::group_by(!!rlang::sym("REP"), !!rlang::sym("TYPE"), 
            	!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
            dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES")))
        reps <- reps %>% 
            dplyr::filter(!!rlang::sym("TO") != "N") %>%
            dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("REP")) %>% 
            dplyr::mutate(PROP = !!rlang::sym("SWITCHES")/
            	sum(!!rlang::sym("SWITCHES"))) %>%
            dplyr::select(-!!rlang::sym("SWITCHES")) %>% 
            tidyr::spread(!!rlang::sym("TYPE"), !!rlang::sym("PROP"))
    }else{        
        reps <- switches  %>%
            dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("CLONE"), 
            	!!rlang::sym("REP")) %>% 
            dplyr::mutate(PROP = !!rlang::sym("SWITCHES")/
            	sum(!!rlang::sym("SWITCHES"))) %>% 
            dplyr::select(-!!rlang::sym("SWITCHES")) %>% 
            tidyr::spread(!!rlang::sym("TYPE"), !!rlang::sym("PROP"))
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
            !reps$TO %in% remove, ]
    }

    if(!bylineage){
        if(alternative[1] == "two.sided"){
            means <- reps %>% 
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>%
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "greater"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),
                    DELTA = mean(!!rlang::sym("DELTA")))
        }            
    }else{
        if(alternative[1] == "two.sided"){
            means <- reps %>% 
                dplyr::group_by(!!rlang::sym("CLONE"), 
                    !!rlang::sym("FROM"), !!rlang::sym("TO")) %>%
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "greater"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("CLONE"),
                    !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("CLONE"),
                    !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),
                    DELTA = mean(!!rlang::sym("DELTA")))
        }                            
    }

    means$STAT <- "SP"
    reps$STAT <- "SP"
    means$REPS <- length(unique(reps$REP))
    results <- list()
    results$reps <- reps
    results$means <- means
    results
}

#' Performs SC (switch count) test on switch data
#' 
#' \code{testSC} performs an SC test
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
#' @seealso Uses output from \link{bootstrapTrees}. Related to \link{testPS}
#' and \link{testSP}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleDb)
#' ExampleDb$sample_id = sample(ExampleDb$sample_id)
#' clones = formatClones(ExampleDb, trait="sample_id")
#' btrees = bootstrapTrees(clones[1:2], bootstraps=100, nproc=1,
#'    igphyml=igphyml, trait="sample_id", id="temp", dir="temp")
#' testSC(btrees$switches)
#' }
#' @export
testSC <- function(switches,dropzeros=TRUE,
	bylineage=FALSE, pseudocount=0, from=NULL, to=NULL,
	permuteAll=FALSE, alternative=c("two.sided","greater","less")){

    permute <-dplyr::quo(!!rlang::sym("PERMUTE"))
    if(permuteAll){
        permute <-dplyr::quo(!!rlang::sym("PERMUTEALL"))
    }

    switches <- switches %>% 
        dplyr::filter(!!rlang::sym("TO") != !!rlang::sym("FROM") & 
        	!!rlang::sym("FROM") != "UCA")

    if(!is.null(from)){
        from <-dplyr::enquo(from)
        switches <- dplyr::filter(switches, !!rlang::sym("FROM") == !!from)
    }
    if(!is.null(to)){
        to <-dplyr::enquo(to)
        switches <- dplyr::filter(switches, !!rlang::sym("TO") == !!to)
    }

    if(!bylineage){
        reps <- switches  %>%
            dplyr::group_by(!!rlang::sym("REP"), !!rlang::sym("TYPE"), 
            	!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
            dplyr::summarize(SWITCHES = sum(!!rlang::sym("SWITCHES")))

        reps <- reps %>% 
            dplyr::filter(!!rlang::sym("TO") != "N") %>%
            dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("REP")) %>% 
            dplyr::mutate(COUNT = !!rlang::sym("SWITCHES")) %>%
            dplyr::select(-!!rlang::sym("SWITCHES")) %>% 
            tidyr::spread(!!rlang::sym("TYPE"), !!rlang::sym("COUNT"))
    }else{
        reps <- switches  %>%
            dplyr::group_by(!!rlang::sym("TYPE"), !!rlang::sym("CLONE"), 
            	!!rlang::sym("REP")) %>% 
            dplyr::mutate(COUNT = !!rlang::sym("SWITCHES")) %>% 
            dplyr::select(-!!rlang::sym("SWITCHES")) %>% 
            tidyr::spread(!!rlang::sym("TYPE"), !!rlang::sym("COUNT"))
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
            !reps$TO %in% remove, ]
    }

    if(!bylineage){
        if(alternative[1] == "two.sided"){
            means <- reps %>% 
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>%
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "greater"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),
                    DELTA = mean(!!rlang::sym("DELTA")))
        }            
    }else{
        if(alternative[1] == "two.sided"){
            means <- reps %>% 
                dplyr::group_by(!!rlang::sym("CLONE"), 
                    !!rlang::sym("FROM"), !!rlang::sym("TO")) %>%
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0)*0.5 + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "greater"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("CLONE"),
                    !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("CLONE"),
                    !!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!rlang::sym("PERMUTE")),
                    PLT = (sum(!!rlang::sym("DELTA") > 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),
                    DELTA = mean(!!rlang::sym("DELTA")))
        }                            
    }
    means$STAT <- "SC"
    reps$STAT <- "SC"
    means$REPS <- length(unique(reps$REP))
    results <- list()
    results$reps <- reps
    results$means <- means
    results
}

#' Performs root-to-tip regression test on set of trees
#' 
#' \code{rootToTop} performs root-to-tip regression permutation test
#' @param    trees        Tibble with trees and data 
#' @param    time         Column id for timepoint
#' @param    permutations Number of permutations for test
#' @param    germline     Germline sequence name
#' @param    minlength    Branch lengths to collapse in trees
#' @param    alternative   Perform test for each lineage individually?
#'
#' @return   A \code{tibble} with pearson correlation between divergene
#' and time, mean permuted correlation, p value(s), number of permutations,
#' and number of sequences
#'
#' @details
#' Output data table columns:
#' clone_id = clone id
#' observed = observed pearson correlation
#' permuted = mean permuted correlation
#' pgt = p value for DELTA < 0
#' plt = p value for DELTA > 0
#'  
#' @seealso Uses output from \link{getTrees}.
#' @export
rootToTip <- function(trees, time="time", permutations=1000,
    germline="Germline", minlength=0.001,
    alternative=c("two.sided","greater","less")){

    # perform root-tip regressions
    regressions <- tibble()
    for(cloneid in unique(trees$clone_id)){
        print(cloneid)
        temp <- dplyr::filter(trees,clone_id == cloneid)
        tree <- temp$trees[[1]]
        data <- temp$data[[1]]@data
    
        if(n_distinct(data[[time]]) == 1 || 
            n_distinct(data[[time]]) == 1){
            next
        }

        dseq = data[is.na(data[[time]]),]$sequence_id
        if(length(dseq) > 0){
            data <- dplyr::filter(data,!sequence_id %in% dseq)
            tree <- ape::drop.tip(tree,tip=dseq)
        }
        
        tips <- tree$tip.label
        tree$edge.length[tree$edge.length < minlength] <- 0
        tree <- ape::di2multi(tree,tol=minlength)
        uca <- ape::getMRCA(tree,tip=tips[!grepl(germline,tips)])
        co <- ape::dist.nodes(tree)
        dist <- co[1:length(tree$tip.label),uca]
        names(dist) <- tree$tip.label
        dist <- dist[!grepl(germline,names(dist))]
    
        # add cophenetic distance to data tibble
        data$divergence <- dist[data$sequence_id]
    
        # get observed and permuted correlation between divergence and time
        observed_cor <- cor(data$divergence,data[[time]])
        perm_temp <- data
        perm_cor <- rep(1,length=permutations)
        for(p in 1:permutations){
            perm_temp[[time]] <- sample(data[[time]],replace=FALSE)
            perm_cor[p] <- cor(perm_temp$divergence,perm_temp[[time]])
        }
    
        # collect results
        results <- tibble(
            clone_id=cloneid,
            observed=observed_cor,
            permuted=mean(perm_cor),
            pv_gt = sum(perm_cor >= observed_cor)/permutations,
            pv_lt = sum(perm_cor <= observed_cor)/permutations,
            nperm = permutations,
            nseq = nrow(data)
            )
        regressions <- bind_rows(regressions,results)
    }
    return(regressions)
}