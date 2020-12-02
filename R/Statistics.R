# Summary statistics and data transformation for switch count distributions

#' Performs PS (parsimony score) test on switch data
#' 
#' \code{testPS} performs a PS test
#' @param    switches     Data frame from bootstrapTrees
#' @param    bylineage    Perform test for each lineage individually? (FALSE)
#' @param    pseudocount  Pseudocount for P value calculations
#' @param    alternative  Perform one-sided (\code{greater} or \code{less})
#'                          or \code{two.sided} test
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
#' @param    alternative  Perform one-sided (\code{greater} or \code{less})
#'                          or \code{two.sided} test
#' @param    binom        Calculate binomial p value
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
    bylineage=FALSE, pseudocount=0, alternative=c("two.sided","greater","less"),
    binom=FALSE){

    if(binom){
        if(!bylineage){
            warning("binom=TRUE, setting bylineage to TRUE")
            bylineage = TRUE
        }
        if(alternative[1] != "greater"){
            warning("binom=TRUE, setting alternative to greater")
            alternative = "greater"
        }
    }

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
            dplyr::filter(!!rlang::sym("TO") != "N") %>%
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
            dplyr::pull(rlang::sym("FROM"))
        
        to_type <- 
            reps %>% dplyr::group_by(!!rlang::sym("TO")) %>%
            dplyr::summarize(d = sum(!!rlang::sym("RECON"))) %>%
            dplyr::filter(!!rlang::sym("d") == 0) %>%
            dplyr::pull(rlang::sym("TO"))
        
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

    if(binom){
        if(dropzeros){
            means = means %>%
                filter(!!rlang::sym("RECON") != 0 | !!rlang::sym("PERMUTE") != 0)
        }
        means = means %>%
            dplyr::group_by(!!rlang::sym("FROM"),!!rlang::sym("TO")) %>%
            summarize(CLONES=n(),
                POSITIVE=sum(!!rlang::sym("DELTA") > 0),
                P=stats::binom.test(!!rlang::sym("POSITIVE"),
                    !!rlang::sym("CLONES"),
                    alternative="greater")$p.value)
        means$TEST = "BINOM"
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
#' @param    alternative  Perform one-sided (\code{greater} or \code{less})
#'                          or \code{two.sided} test
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
            dplyr::filter(!!rlang::sym("TO") != "N") %>%
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
            dplyr::pull(rlang::sym("FROM"))
        
        to_type <- 
            reps %>% dplyr::group_by(!!rlang::sym("TO")) %>%
            dplyr::summarize(d = sum(!!rlang::sym("RECON"))) %>%
            dplyr::filter(!!rlang::sym("d") == 0) %>%
            dplyr::pull(rlang::sym("TO"))
        
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


#' Get divergence from MRCA for each tip
#' 
#' \code{rootToTop} performs root-to-tip regression permutation test
#' @param    phy          Tree object
#' @param    minlength    Branch lengths to collapse in trees
#' @return   A named vector of each tip's divergence from the tree's MRCA.
#'
#' @export
getDivergence = function(phy, minlength=0.001){
    tips <- phy$tip.label
    phy$edge.length[phy$edge.length < minlength] <- 0
    phy <- ape::di2multi(phy,tol=minlength)
    uca <- ape::getMRCA(phy,tip=tips)
    co <- ape::dist.nodes(phy)
    dist <- co[1:length(tips),uca]
    names(dist) <- tips
    dist
}

#' Resolve polytomies to have the minimum number of single timepoint clades
#' 
#' @param    phy          Tree object
#' @param    clone        airrClone data object corresponding to \code{phy}
#' @param    time         Column name holding numeric time information
#' @param    sequence     Column name holding sequence ID
#' @param    germline     Germline sequence name
#' @param    minlength    Branch lengths to collapse in trees
#' @param    verbose      Print lots of rubbish while running?
#' @return   A \code{phylo} tree object in which polytomies are resolved to 
#' have the minimum number of single timepoint clades.
#'
#' @details
#' Iteratively identifies polytomies (clusters of < minlength branches),
#' prunes each descendant branch, combines clades with the same timepoint
#' before grouping them back together. Checks to make sure that the divergence
#' of each tip is the same after resolution.
#'  
#' @seealso Uses output from \link{getTrees} during \link{correlationTest}.
#' @export
resolvePolytomies = function(phy, clone, minlength=0.001,
    time="time", sequence="sequence_id", germline = "Germline",
    verbose=FALSE){
    
    data <- clone@data
    phy <- rerootTree(di2multi(phy,tol=minlength),
        germline=germline)
    tips <- phy$tip.label
    if(sum(!data[[sequence]] %in% phy$tip.label) > 0){
        stop("Tree and data sequence ids don't match")
    }

    ittybitty_branch = 1e-7
    # ape does weird things with zero branches
    phy$edge.length[phy$edge.length == 0] = ittybitty_branch

    odivergence <- getDivergence(phy,minlength)
    polytomies <- table(phy$edge[,1])
    polytomies <- names(polytomies[polytomies > 2])
    while(length(polytomies) > 0){
        if(verbose){
            print(polytomies)
        }
        # for first polytomy in list, remove each subclade
        # and record the time of its sequences
        target_node <- as.numeric(polytomies[1])
        node_index <- which(phy$edge[,1] == target_node)
        node_ns <- phy$edge[node_index,2]
        edge_ls <- phy$edge.length[node_index]
        clades <- list()
        times <- list()
        for(i in 1:length(node_ns)){
            nc <- node_ns[i]
            # unfortuantely processing is different if clade is > 1 tip.
            if(nc <= length(phy$tip.label)){
                clades[[i]] <- ape::keep.tip(phy, tip = nc)
                times[[i]] <- unique(data[data[[sequence]] %in%
                    clades[[i]]$tip.label,][[time]])
                clades[[i]]$root.edge <- edge_ls[[i]] #store root edge
                clades[[i]]$edge.length <- 0
            }else{
                tips <- ape::extract.clade(phy, node = nc, 
                    collapse.singles=TRUE)$tip.label
                clades[[i]] <- ape::keep.tip(phy, tip = tips)
                clades[[i]]$root.edge <- edge_ls[[i]] #store root edge
                times[[i]] <- unique(data[data[[sequence]] %in%
                    clades[[i]]$tip.label,][[time]])
                if(length(times[[i]]) > 1){
                    times[[i]] <- "multi"
                }
            }
        }

        # combine all clades of each time type in a balanced
        # manner, so they remain clades
        topclades <- list()
        distinct_times <- unique(unlist(times))
        for(ctime in distinct_times){
            tclades <- clades[unlist(lapply(times,
                function(x)x==ctime))]
            clade <- tclades[[1]]
            if(length(tclades) > 1){
                for(i in 2:length(tclades)){
                    clade <- ape::bind.tree(clade,tclades[[i]],
                        position=clade$root.edge)
                    clade <- ape::collapse.singles(clade)
                    clade$root.edge <- ittybitty_branch
                }
            }
            topclades[[as.character(ctime)]] <- clade
        }

        # combine top clades into a polytomy
        polytomy <- topclades[[1]]
        if(length(topclades) > 1){
            for(i in 2:length(topclades)){
                polytomy <- ape::bind.tree(polytomy,topclades[[i]],
                            position=polytomy$root.edge)
                polytomy <- ape::collapse.singles(polytomy)
                polytomy$root.edge <- ittybitty_branch
            }
        }
        tphy <- ape::drop.tip(phy, collapse.singles=TRUE, trim.internal=FALSE, tip = 
            ape::extract.clade(phy, node = target_node)$tip.label, subtree=TRUE)
        dumb_tip <- paste0("[",length(polytomy$tip.label),"_tips]")
        ntree <- ape::bind.tree(tphy, polytomy, where=which(tphy$tip.label==dumb_tip),
            position=polytomy$root.edge)
        ntree <- ape::drop.tip(ntree, tip=dumb_tip)
        phy <- ntree
        polytomies <- table(phy$edge[,1])
        polytomies <- names(polytomies[polytomies > 2])
    }

    ndivergence <- getDivergence(phy,minlength)

    if(length(ndivergence) != length(odivergence)){
        stop("Number of tips corrupted during polytomy resolution")
    }
    if(sum(odivergence[names(ndivergence)] != ndivergence)){
        stop("Divergence corrupted during polytomy resolution")   
    }
    phy
}

#' Resolve polytomies to have the minimum number of single timepoint clades
#' 
#' \code{rootToTop} performs root-to-tip regression permutation test
#' @param    phy          Tree object
#' @param    clone        airrClone data object corresponding to \code{phy}
#' @param    permutations Number of permutations to run
#' @param    polyresolve  Resolve polytomies to have a minimum number of 
#'                         single timepoint clades
#' @param    permutation  Permute among single timepoint clades or uniformly
#'                         among tips
#' @param    time         Column name holding numeric time information
#' @param    sequence     Column name holding sequence ID
#' @param    germline     Germline sequence name
#' @param    minlength    Branch lengths to collapse in trees
#' @param    verbose      Print lots of rubbish while running?
#' @param    alternative  Is alternative that the randomized correlation are greater than 
#'                         or equal to observed, or greater/less than?
#' @return   A list of statistics from running the permutation test.
#'
#' @details
#'  See \link{correlationTest} for details
#' @seealso \link{correlationTest}.
runCorrelationTest = function(phy, clone, permutations, minlength=0.001,
    polyresolve = TRUE, permutation = c("clustered", "uniform"), 
    time="time", sequence="sequence_id", germline = "Germline",
    verbose=TRUE, alternative = c("greater","two.sided")){

    data <- clone@data
    if(verbose){
        print(paste("Analyzing clone: ",clone@clone))
    }

    if(polyresolve){
        if(verbose){
            print("resolving polytomies to single timepoint clades")
        }
        phy <- resolvePolytomies(phy, clone, minlength,
            time=time, sequence=sequence, germline = germline,
            verbose=verbose)
    }

    phy <- ape::drop.tip(phy,germline)
    tips <- phy$tip.label
    if(sum(!data[[sequence]] %in% phy$tip.label) > 0){
        stop("Tree and data sequence ids don't match")
    }

    # dates is named and ordered by tree tips
    dates <- data[[time]]
    names(dates) <- data[[sequence]]
    dates <- dates[tips]

    divergence <- getDivergence(phy, minlength)
    data$divergence <- divergence[data[[sequence]]]

    cl <- 1:length(dates) # each tip is its own cluster
    if(permutation[1] == "clustered"){
        # define all monophyletic clades with > 1 tips
        gc <- lapply(ape::subtrees(phy),function(x)x$tip.label)
        # identify which clades are single timepoint
        single.date.clades <- which(unlist(lapply(gc,function(x)
            dplyr::n_distinct(dates[x])))==1)
        if(length(single.date.clades)>0){
            for(i in 1:length(single.date.clades)){
                m <- match(gc[[single.date.clades[i]]],phy$tip.label)
                cl[m] <- cl[m[1]] # assign clusters to be the same across single timepoint clades
            }
        }
        cl <- match(cl,unique(cl)) #make clusters 1:(number of clusters)
        if(length(unique(cl))==1) stop("Only one cluster in data.")
    }
    # assign clusters to each sequence
    names(cl) <- tips
    m <- match(data[[sequence]],names(cl))
    data$cluster <- cl[data$sequence_id]

    counts <- table(data$cluster,data$time)
    if(sum(rowSums(counts > 0) > 1) > 0){
        stop("Clusters are not single timepoint!")
    }

    # get ordered list of times for each cluster
    ctimes <- unlist(lapply(sort(unique(data$cluster)), function(x)
        unique(data[data$cluster == x,][[time]])))

    true <- stats::cor(data$time,data$divergence)
    random <- rep(0, length=permutations)
    for(i in 1:permutations){
        times <- ctimes[sample(1:length(ctimes))]
        data$random_time <- times[data$cluster]
        random[i] <- stats::cor(data$divergence, data$random_time)
    }

    nclust <- dplyr::n_distinct(cl)
    nposs <-exp(lfactorial(nclust)-sum(lfactorial(table(ctimes))))

    minp <- 1
    if(alternative[1] == "two.sided"){
        minp <- 0.5
    }

    clone@data <- data
    slope <- summary(stats::lm(data$divergence ~ data$time))$coefficients[2,1]
    results <- list(correlation=true)
    results[["clone"]] <- clone
    results[["tree"]] <- phy
    results[["random"]] <- random
    results[["random_correlation"]] <- mean(random)
    random <- c(random,true)
    results[["p_gt"]] <- (sum(true < random) + 
        sum(true == random)*0.5)/length(random)
    results[["p_lt"]] <- (sum(true > random) + 
        sum(true == random)*0.5)/length(random)
    results[["nposs"]] <- nposs
    results[["nclust"]] <- nclust
    results[["p"]] <- mean(true <= random)
    results[["min_p"]] <- max(minp/nposs, minp/(permutations+1))
    results[["slope"]] <- slope
    results
}

#' Run date randomization test for temporal signal on a set of trees.
#' 
#' \code{correlationTest} performs root-to-tip regression date randomization test
#' @param    clones       A \code{tibble} object containing airrClone and \code{phylo} objects
#' @param    permutations Number of permutations to run
#' @param    polyresolve  Resolve polytomies to have a minimum number of 
#'                         single timepoint clades
#' @param    permutation  Permute among single timepoint clades or uniformly
#'                         among tips
#' @param    time         Column name holding numeric time information
#' @param    sequence     Column name holding sequence ID
#' @param    germline     Germline sequence name
#' @param    minlength    Branch lengths to collapse in trees
#' @param    verbose      Print lots of rubbish while running?
#' @param    storeTree    Store the tree used?
#' @param    alternative  Is alternative that the randomized correlation are greater than 
#'                         or equal to observed, or greater/less than?
#' @param    nproc        Number of cores to use for calculations. Parallelizes by tree.
#' @return   A \code{tibble} with the same columns as clones, but additional
#' columns corresponding to test statistics for each clone. 
#'
#' @details
#'  Object returned contains these columns which are added or modified from input:
#'  \itemize{
#'   \item  \code{data}: airrClone object, same as input but with additional columns 
#'        "cluster" which correspond to permutation cluster, and "divergence."
#'   \item  \code{slope}: Slope of linear regression between divergence and time.
#'   \item  \code{correlation}: Correlation between divergence and time.
#'   \item  \code{p}: p value of correlation compared to permuted correlations.
#'   \item  \code{random_correlation}: Mean correlation of permutation replicates.
#'   \item  \code{min_p}: Minimum p value of data, determined by either the number of
#'         distinct clade/timepoint combinations or number of permutations.
#'   \item  \code{nposs}: Number of possible distinct timepoint/clade combinations.
#'   \item  \code{nclust}: Number of clusters used in permutation. If permutation="uniform"
#'         this is the number of tips.
#'   \item  \code{p_gt/p_lt}: P value that permuted correlations are greater or less 
#'         than observed correlation. Only returned if alternative = "two.sided"
#'   \item  \code{test_trees}:  The \link{phylo} tree objects used, possibly with resolved polytomies.
#' }
#' @seealso Uses output from \code{getTrees}.
#' @export
correlationTest = function(clones, permutations=1000, minlength=0.001,
    permutation = c("clustered", "uniform"), time="time", 
    sequence="sequence_id", germline = "Germline",
    verbose=FALSE, polyresolve = TRUE,
    alternative = c("greater","two.sided"),
    storeTree = FALSE, nproc=1){

    results <- parallel::mclapply(1:nrow(clones),function(x)
        runCorrelationTest(clones$trees[[x]], clones$data[[x]],
        permutations, minlength=minlength, polyresolve=polyresolve,
        permutation=permutation, time= time, 
        sequence=sequence, germline=germline,
        verbose=verbose, alternative=alternative),
        mc.cores=nproc)

    clones$data <- lapply(results,function(x)x$clone)
    if(storeTree){
        clones$test_trees <- lapply(results,function(x)x$tree)
    }

    cols <- c("slope", "p", "correlation", "random_correlation",
        "min_p", "nposs", "nclust")
    if(alternative[1] == "two.sided"){
        cols <- c(cols, c("p_gt", "p_lt"))
    }
    for(col in cols){
        clones[[col]] <- unlist(lapply(results,function(x)x[[col]]))
    }
    clones
}