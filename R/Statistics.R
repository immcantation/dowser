# Summary statistics and data transformation for switch count distributions

#' Performs PS (parsimony score) test on switch data
#' 
#' \code{testPS} performs a PS test
#' @param    switches     Data frame from findSwitches
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
#' @seealso Uses output from \link{findSwitches}. Related to \link{testSP}
#' and \link{testSC}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleAirr)
#' ExampleAirr$sample_id <- sample(ExampleAirr$sample_id)
#' clones <- formatClones(ExampleAirr, trait="sample_id")
#' btrees <- findSwitches(clones[1:2], bootstraps=10, nproc=1,
#'    igphyml=igphyml, trait="sample_id")
#' testPS(btrees$switches)
#' }
#' @export
testPS <- function(switches, bylineage=FALSE, pseudocount=0,
    alternative=c("less","two.sided","greater")){
    switches <- switches %>% 
        dplyr::filter(!!rlang::sym("TO") != "N" & 
            !!rlang::sym("TO") != !!rlang::sym("FROM") &
             !!rlang::sym("FROM") != "UCA")
    if(nrow(switches) == 0){
        stop("No switches left after filtering!")
    }
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
#' @param    switches     Data frame from findSwitches
#' @param    permuteAll   Permute among trees?
#' @param    from         Include only switches from this state?
#' @param    to           Include only switches to this state?
#' @param    dropzeroes    Drop switches with zero counts?
#' @param    bylineage    Perform test for each lineage individually?
#' @param    pseudocount  Pseudocount for P value calculations
#' @param    alternative  Perform one-sided (\code{greater} or \code{less})
#'                          or \code{two.sided} test
#' @param    tip_switch   maximum tip/switch ratio
#' @param    exclude      exclude clones with tip/switch ratio > \code{tip_switch}?
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
#' @seealso Uses output from \link{findSwitches}. Related to \link{testPS}
#' and \link{testSC}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleAirr)
#' ExampleAirr$sample_id = sample(ExampleAirr$sample_id)
#' clones = formatClones(ExampleAirr, trait="sample_id")
#' btrees = findSwitches(clones[1:2], bootstraps=10, nproc=1,
#'    igphyml=igphyml, trait="sample_id")
#' testSP(btrees$switches)
#' }
#' @export
testSP <- function(switches, permuteAll=FALSE, 
    from=NULL, to=NULL, dropzeroes=TRUE,
    bylineage=FALSE, pseudocount=0, alternative=c("greater","two.sided","less"),
    tip_switch=20, exclude=FALSE){

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

    switches <- switches %>% 
        dplyr::filter(!!rlang::sym("TO") != "N")

    if(nrow(switches) == 0){
        stop("No switches left after filtering!")
    }

    if(exclude){
        # separate tips
        tips <- dplyr::filter(switches, !!rlang::sym("TO")=="NTIP" &
            !!rlang::sym("TYPE")=="RECON")
    
        # germline doesn't count in downsampling algorithm
        tips$SWITCHES <- tips$SWITCHES - 1
    
        counts <- testPS(switches, bylineage=TRUE)$means
        m <- match(counts$CLONE, tips$CLONE)
        counts$tips <- tips[m,]$SWITCHES
        counts$ratio <- counts$tips/counts$RECON
    
        excluded <- dplyr::filter(counts, !!rlang::sym("ratio") > tip_switch)$CLONE
        if(length(excluded) > 0){
            warning(paste("Excluding clone(s)",paste(excluded,collapse=",")
                ,"due to high tip/switch ratio"))
    
            switches <- dplyr::filter(switches,!(!!rlang::sym("CLONE") %in% excluded))
        }
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
    
    if(dropzeroes){
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

    if(sum(is.na(reps$DELTA)) > 0){
        na_reps <- unique(reps[is.na(reps$DELTA),]$REP)
        warning(paste0("NA delta values in ",
            length(na_reps)," replicates, discarding"))
        reps <- dplyr::filter(reps, !(!!rlang::sym("REP") %in% na_reps))
    }

    if(!bylineage){
        if(alternative[1] == "two.sided"){
            means <- reps %>% 
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>%
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
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
#' @param    switches     Data frame from findSwitches
#' @param    permuteAll   Permute among trees?
#' @param    from         Include only switches from this state?
#' @param    to           Include only switches to this state?
#' @param    dropzeroes    Drop switches with zero counts?
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
#' @seealso Uses output from \link{findSwitches}. Related to \link{testPS}
#' and \link{testSP}.
#' @examples
#' \dontrun{
#' igphyml <- "~/apps/igphyml/src/igphyml"
#' data(ExampleAirr)
#' ExampleAirr$sample_id = sample(ExampleAirr$sample_id)
#' clones = formatClones(ExampleAirr, trait="sample_id")
#' btrees = findSwitches(clones[1:2], bootstraps=100, nproc=1,
#'    igphyml=igphyml, trait="sample_id", id="temp", dir="temp")
#' testSC(btrees$switches)
#' }
#' @export
testSC <- function(switches,dropzeroes=TRUE,
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

    switches <- switches %>% 
        dplyr::filter(!!rlang::sym("TO") != "N")

    if(nrow(switches) == 0){
        stop("No switches left after filtering!")
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

    if(dropzeroes){
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
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
                    PGT = (sum(!!rlang::sym("DELTA") < 0) + 
                        sum(!!rlang::sym("DELTA") == 0) + pseudocount)/
                        (dplyr::n() + pseudocount),    
                    DELTA = mean(!!rlang::sym("DELTA")))
        }else if(alternative[1] == "less"){
            means <- reps %>%
                dplyr::group_by(!!rlang::sym("FROM"), !!rlang::sym("TO")) %>% 
                dplyr::summarize(
                    RECON = mean(!!rlang::sym("RECON")),
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
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
                    PERMUTE = mean(!!permute),
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


#' Get divergence from root of tree for each tip
#' 
#' \code{getDivergence} get sum of branch lengths leading from the 
#' root of the tree. If the germline sequence is included in the tree,
#' this will equal the germline divergence. If germline removed,
#' this will equal the MRCA divergence
#' @param    phy          Tree object
#' @param    minlength    Branch lengths to collapse in trees
#' @return   A named vector of each tip's divergence from the tree's root.
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
    phy <- rerootTree(ape::di2multi(phy,tol=minlength),
        germline=germline,verbose=0)
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
            if(nc <= length(phy$tip.label)){ #clade is a tip
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
        pruned_tip <- paste0("[",length(polytomy$tip.label),"_tips]")
        ntree <- ape::bind.tree(tphy, polytomy, where=which(tphy$tip.label==pruned_tip),
            position=polytomy$root.edge)
        ntree <- ape::drop.tip(ntree, tip=pruned_tip)
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

#' Run correlationTest, based on https://doi.org/10.1111/2041-210X.12466
#' 
#' \code{runCorrelationTest} performs root-to-tip regression permutation test
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

    if(dplyr::n_distinct(dates)==1){
        warning(paste("Only one timepoint cluster in clone",clone@clone))
        results <- list(correlation=NA)
        results[["clone"]] <- clone
        results[["tree"]] <- phy
        results[["random"]] <- NA
        results[["random_correlation"]] <- NA
        results[["p_gt"]] <- NA
        results[["p_lt"]] <- NA
        results[["nposs"]] <- NA
        results[["nclust"]] <- NA
        results[["p"]] <- NA
        results[["min_p"]] <- NA
        results[["slope"]] <- NA
        return(results)
    }

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
    }
    # assign clusters to each sequence
    names(cl) <- tips
    m <- match(data[[sequence]],names(cl))
    data$cluster <- cl[data$sequence_id]

    counts <- table(data$cluster,data[[time]])
    if(sum(rowSums(counts > 0) > 1) > 0){
        stop("Clusters are not single timepoint!")
    }

    # get ordered list of times for each cluster
    ctimes <- unlist(lapply(sort(unique(data$cluster)), function(x)
        unique(data[data$cluster == x,][[time]])))

    true <- stats::cor(data[[time]],data$divergence)
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
    slope <- summary(stats::lm(data$divergence ~ data[[time]]))$coefficients[2,1]
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
#'                        single timepoint clades
#' @param    perm_type    Permute among single timepoint clades or uniformly
#'                        among tips
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
#'   \item  \code{nclust}: Number of clusters used in permutation. If perm_type="uniform"
#'         this is the number of tips.
#'   \item  \code{p_gt/p_lt}: P value that permuted correlations are greater or less 
#'         than observed correlation. Only returned if alternative = "two.sided"
#'   \item  \code{test_trees}:  The \link{phylo} tree objects used, possibly with
#'          resolved polytomies.
#' }
#' @seealso Uses output from \code{getTrees}.
#' @export
correlationTest = function(clones, permutations=1000, minlength=0.001,
    perm_type = c("clustered", "uniform"), time="time", 
    sequence="sequence_id", germline = "Germline",
    verbose=FALSE, polyresolve = TRUE,
    alternative = c("greater","two.sided"),
    storeTree = FALSE, nproc=1){

    if(!"tbl" %in% class(clones)){
        print(paste("clones is of class",class(clones)))
        stop("clones must be a tibble of airrClone objects!")
    }else{
        if(!inherits(clones$data[[1]], "airrClone")){
            print(paste("clones is list of class",class(clones$data[[1]])))
            stop("clones must be a list of airrClone objects!")
        }
    }
    if(!"trees" %in% names(clones)){
        stop("clones must have trees column!")
    }
    time_check <- unlist(lapply(clones$data, function(x)!time %in% names(x@data)))
    if(sum(time_check) > 0){
        stop(paste("Time column",time,"not found in clone object (must be trait value in formatClones)"))
    }
    time_check <- unlist(lapply(clones$data, function(x)!is.numeric(x@data[[time]])))
    if(sum(time_check) > 0){
        stop(paste("Time column",time,"contains non-numeric values, impossible to continue"))
    }

    results <- parallel::mclapply(1:nrow(clones),function(x)
        runCorrelationTest(clones$trees[[x]], clones$data[[x]],
        permutations, minlength=minlength, polyresolve=polyresolve,
        permutation=perm_type, time= time, 
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

#' Finds the Robinson-Fould's cluster distance between phylogenies. 
#' 
#' \code{calcRF} Calculates the RF distance between two phylogenetic trees with 
#'               the same tips and tip labels.
#' @param tree_1         A \code{phylo} object
#' @param tree_2         A \code{phylo} object
#' @param nproc          Number of cores to use for calculations.
#'
#' @return   The RF cluster value for the two input trees.
#'  
#' @export
calcRF <- function(tree_1, tree_2, nproc = 1){
  tip_amount_check <- length(tree_1$tip.label) == length(tree_2$tip.label)
  if(!tip_amount_check){
    stop("trees do not have the same number of tips")
  }
  tip_check <- dplyr::setdiff(tree_1$tip.label, tree_2$tip.label)
  # change this 
  if(!identical(tip_check, character(0))){
    stop("tree tip labels are not identical")
  }
  tree_1_df <- splits_func(list(tree_1),1)
  tree_2_df <- splits_func(list(tree_2), 1)
  
  total_mismatches <- unlist(parallel::mclapply(1:nrow(tree_1_df), function(x){
    tree_1_sub <- tree_1_df$found[[x]]
    mismatch_vector <- unlist(lapply(1:nrow(tree_2_df), function(y){
      mismatches_1 <- dplyr::setdiff(tree_2_df$found[[y]], tree_1_sub)
      mismatches_2 <- dplyr::setdiff(tree_1_sub, tree_2_df$found[[y]])
      if(identical(mismatches_1, character(0)) & identical(mismatches_2, character(0))){
        value <- "match"      
      } else{
        value <- "mismatch"
      }
      return(value)
    }))
    if("match" %in% mismatch_vector){
      tobind <- 0
    } else{
      tobind <- 1
    }
    return(tobind)
  }, mc.cores = nproc))
  total_mismatches <- sum(total_mismatches)

  clone_mismatches <- unlist(parallel::mclapply(1:nrow(tree_2_df), function(x){
    tree_2_sub <- tree_2_df$found[[x]]
    mismatch_vector <- unlist(lapply(1:nrow(tree_1_df), function(y){
      mismatches_1 <- setdiff(tree_1_df$found[[y]], tree_2_sub)
      mismatches_2 <- setdiff(tree_2_sub, tree_1_df$found[[y]])
      if(identical(mismatches_1, character(0)) & identical(mismatches_2, character(0))){
        value <- "match"
      } else{
        value <- "mismatch"
      }
      return(value)
    }))
    if("match" %in% mismatch_vector){
      tobind <- 0
    } else{
      tobind <- 1
    }
    return(tobind)
  }, mc.cores = nproc))
  clone_mismatches <- sum(clone_mismatches)
  
  all_mismatches <- total_mismatches + clone_mismatches
  
  return(all_mismatches)
}

#' Get parameter estimates from posterior distributions
#' 
#' \code{getParams} performs root-to-tip regression date randomization test
#' @param  clones    A \code{tibble} object containing \code{airrClone} and 
#'                   parameters_posterior column
#' @param  burnin    percent of initial samples to discard (1-100)
#' @param  tracefile file name for trace plots
#' @param  width     width of plot in inches
#' @param  height    height of plot in inches
#' @param  ...       optional arguments passed to grDevices::pdf
#' @return   A \code{tibble} with the same columns as clones, but additional
#' column parameters corresponding to parameter estimates
#'
#' @export
getParams = function(clones, burnin=10, tracefile=NULL, width=8.5, height=11, ...){

    if(!"tbl" %in% class(clones)){
        print(paste("clones is of class",class(clones)))
        stop("clones must be a tibble of airrClone objects!")
    }else{
        if(!inherits(clones$data[[1]], "airrClone")){
            print(paste("clones is list of class",class(clones$data[[1]])))
            stop("clones must be a list of airrClone objects!")
        }
    }
    if(!"parameters_posterior" %in% names(clones$trees[[1]])){
        stop("tree objects must have parameters_posterior (from BEAST)!")
    }
    if("parameters" %in% names(clones)){
        warning("overwriting parameters column")
    }

    # TODO add ESS calculation!
    parameters <- list()
    for(i in 1:nrow(clones)){
        post <- clones$trees[[i]]$parameters_posterior
        if(nrow(post) == 0){
            stop(paste("parameters_posterior empty for clone", 
                clones$clone_id[[i]]))
        }
        # assumes samples are evenly distributed and that all params
        # are included at each sampling event
        sample_range = range(post$Sample)
        burn <- floor(burnin/100 * (sample_range[2] - sample_range[1]))
        post <- filter(post, !!rlang::sym("Sample") >= burn)
        estimates <- post %>%
            group_by(parameter) %>%
            summarize(estimate = mean(value),
                lci = quantile(value, prob=0.025),
                hci = quantile(value, prob=0.975),
                ess = floor(tryCatch(mcmcse::ess(value),error=function(e)NA)),
                n = n())
        parameters[[i]] <- estimates
    }
    clones$parameters <- parameters

    if(!is.null(tracefile)){
        plotTraces(clones, burnin=burnin, file=tracefile, 
            width=width, height=height, ...)
    }

    clones
}

#' get values for Bayesian Skyline plot
#' 
#' \code{makeSkyline} 
#' @param  logfile   Beast log file
#' @param  treesfile BEAST trees file 
#' @param  burnin    file name for trace plots
#' @param  bins      number of bins for plotting
#' @param  youngest  timepoint of the most recently tip sampled (if 0, backward time used)
#' @param  clone_id  name of the clone being analyzed (if desired)
#' @return   Bayesian Skyline values for given clone
#'
#' @export
makeSkyline <- function(logfile, treesfile, burnin, bins=100, youngest=0, 
    clone_id=NULL){
    
    l <- tryCatch(read.csv(logfile, head=TRUE, sep="\t", comment.char="#"),error=function(e)e)
    if("error" %in% class(l)){
        stop(paste("couldn't open",logfile))
    }
    phylos <- tryCatch(ape::read.nexus(treesfile), error=function(e)e)
    if("error" %in% class(phylos)){
        stop(paste("couldn't open", treesfile))
    }
    params <- tidyr::gather(l, "parameter", "value", -Sample)

    if(!"bPopSizes.1" %in% unique(params$parameter)){
        stop(paste("log file doesn't have pop sizes.",
            "Was it run with skyline tree_prior='coalescent_skyline'?"))
    }

    burn <- floor(length(phylos)*burnin/100)
    samples <- unique(params$Sample)
    if(burn > 0){
      phylos <- phylos[(burn+1):length(phylos)]
      samples <- samples[(burn+1):length(samples)]
      params <- filter(params, Sample %in% samples)
    }
    if(n_distinct(params$Sample) != length(phylos)){
      stop("Parameter and tree posteriors not same length")
    }

    groups <- filter(params, grepl("GroupSizes", parameter))
    pops <- filter(params, grepl("PopSizes", parameter))

    if(sum(pops$value < 0) > 0){
        stop(paste(logfile, "found popsizes < 0, can't continue"))
    }
    if(sum(groups$value < 0) > 0){
        stop(paste(logfile, "found groupsizes < 0, can't continue"))
    }

    pops$index <- as.numeric(gsub("bPopSizes\\.","",pops$parameter))
    groups$index <- as.numeric(gsub("bGroupSizes\\.","",groups$parameter))

    # smallest tree height in log file (this is what tracer seems to do)
    maxheight <- min(filter(params, parameter == "Tree.height")$value)
    if(youngest > 0){
        mintime <- youngest - maxheight;
        maxtime <- youngest;
    }else{
        mintime <- 0
        maxtime <- maxheight - youngest
    }

    binwidth <- (maxtime - mintime)/(bins - 1)

    all_intervals <- tibble()
    for(index in 1:length(phylos)){
      tr <- phylos[[index]]
      sample <- samples[index]
      mrca <- ape::getMRCA(tr, tip=tr$tip.label)
      d <- ape::dist.nodes(tr)
      times <- d[mrca,]
      maxheight <- max(times)
      nodes <- maxheight - times[(length(tr$tip.label)+1):length(times)]
      nodes <- nodes[order(nodes, decreasing=FALSE)]


      pop <- filter(pops, !!rlang::sym("Sample") == sample)
      group <- filter(groups, !!rlang::sym("Sample") == sample)

      temp <- nodes
      results <- tibble()
      for(i in 1:nrow(group)){
        groupsize <- group$value[i]
        popsize <- pop$value[i]
        events <- temp[1:(groupsize)]
        results <- bind_rows(results,
          tibble(end=events[length(events)], interval=i,
            events=length(events), popsize=popsize))
        temp <- temp[-1:-(groupsize)]
      }
      results$sample <- sample
      results$index <- index
      all_intervals <- bind_rows(all_intervals, results)
    }

    indistinct <- all_intervals %>%
        group_by(sample) %>%
        summarize(distinct = dplyr::n_distinct(end),
            n = n()) %>%
        filter(distinct < n) %>%
        pull(sample)

    if(length(indistinct) > 0){
        warning(paste(logfile, "Removing",length(indistinct),
            "samples with indistinct intervals. This shouldn't happen."))
        all_intervals <- dplyr::filter(all_intervals, !sample %in% indistinct)
        if(nrow(all_intervals) == 0){
            stop("No intervals left :-(")
        }
    }
 
    skyline <- tidyr::tibble()
    n_sample <- dplyr::n_distinct(all_intervals$sample)
    bin_intervals <- seq(0, length=bins, by=binwidth)
    interval_bins <- tibble()
    for(j in 1:(length(bin_intervals)-1)){
      bin_start <- bin_intervals[j]
      bin_end <- bin_intervals[j+1]

      matches <- all_intervals %>%
        dplyr::group_by(!!rlang::sym("sample")) %>%
        dplyr::filter(!!rlang::sym("end") > bin_start) %>%
        slice_min(!!rlang::sym("end"))

      if(nrow(matches) != n_sample){
        stop("didn't find some indexes")
      }

      matches <- select(matches, !!rlang::sym("end"), 
        !!rlang::sym("popsize"), !!rlang::sym("sample"), !!rlang::sym("index"))
      matches$bin <- bin_start
      skyline <- bind_rows(skyline, matches)
    }

    # similar behavior to 
    # dr.stat.getQuantile in beast-mcmc
    getQuantiles <- function(x, probs=NULL){
      if(is.null(probs)){
        stop("probs can't be null")
      }
      x <- sort(x)
      index <- sapply(probs, function(y)
        ceiling(length(x)*y))
      x <- x[index]
      names(x) <- paste0(probs*100,"%")
      return(x)
    }

    skyplot <- skyline %>%
      dplyr::group_by(!!rlang::sym("bin")) %>%
      dplyr::summarize(
        mean=mean(!!rlang::sym("popsize")), 
        median=getQuantiles(!!rlang::sym("popsize"), 0.5),
        lci=getQuantiles(!!rlang::sym("popsize"), 0.025), 
        uci=getQuantiles(!!rlang::sym("popsize"), 0.975))

    if(!is.null(clone_id)){
        skyplot$clone_id <- clone_id
    }
    if(youngest != 0){
        skyplot$bin <- youngest - skyplot$bin
    }

    return(skyplot)
}

#' Make data frames for Bayesian skyline plots
#' 
#' \code{makeSkylines} 
#' @param  clones    clone tibble
#' @param  dir       directory of BEAST trees file 
#' @param  time      name of time column
#' @param  bins      number of bins for plotting
#' @param  oldest    age of the oldest tip sampled (if forward time desired)
#' @param  verbose   if 1, print name of clones
#' @param  forward   plot in forward or (FALSE) backward time?
#' @param  nproc     processors for parallelization (by clone)
#' @return   Bayesian Skyline values for given clone
#' @details Burnin set from readBEAST or getTrees
#' @export
getSkylines <- function(clones, dir, time, bins=100, verbose=0, forward=TRUE,
    nproc=1){

    treesfiles <- sapply(clones$data, function(x)
        file.path(dir, paste0(x@clone, ".trees")))

    logfiles <- sapply(clones$data, function(x)
        file.path(dir, paste0(x@clone, ".log")))

    burnins <- sapply(clones$trees, function(x)
        x$burnin)

    if(sum(is.null(burnins)) > 0){
        stop("burnin not found in some tree objects")
    }

    if(forward){
        youngest <- sapply(clones$data, function(x)
            max(as.numeric(x@data[[time]])))
    }else{
        youngest <- rep(0, length=nrow(clones))
    }

    skylines <- parallel::mclapply(1:nrow(clones), function(x){
        if(verbose != 0){
            print(paste(clones$clone_id[x], logfiles[x], 
                treesfiles[x], youngest[x], burnins[x]))
        }
        tryCatch(makeSkyline(logfile=logfiles[x], treesfile=treesfiles[x],
            youngest=youngest[x], burnin=burnins[x], bins=bins, 
            clone_id=clones$clone_id[x]), error=function(e)e)
    }, mc.cores=nproc)

    clones$skyline <- skylines
    for(i in 1:nrow(clones)){
        if("error" %in% class(skylines[[i]])){
            print(skylines[[i]])
            warning(paste("Error making skyline clone,",clones$clone_id[i]))
            clones$skyline[[i]] <- NA
        }
    }
    return(clones)
}


