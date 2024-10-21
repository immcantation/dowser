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
      # CGJ 5/20/24 branch trimming naming was weird with ape update
      # set the node labels to null and it should work the same as previous ape version
      # add check to see if we get more than one tip of the naming pattern
      phy$node.label <- NULL
      tphy <- ape::drop.tip(phy, collapse.singles=TRUE, trim.internal=FALSE, tip = 
                              ape::extract.clade(phy, node = target_node)$tip.label, subtree=TRUE)
      pruned_tip <- paste0("[",length(polytomy$tip.label),"_tips]")
      if(sum(tphy$tip.label == pruned_tip) > 1){
        stop("Polytomy binding point identification failed.")
      }
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
#'   \item  \code{test_trees}:  The \code{phylo} tree objects used, possibly with
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

