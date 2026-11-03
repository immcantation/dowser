# Core of getTimeTreesIterate.
#
# Most functions in this file are internal, except restoreIterateBackup and
# readFromLastCheckpoint, which are intended quality of life improvements.
#
# Runs one BEAST job per clone. When a job ends, its worker slot
# goes to the next clone in the queue at once. No clone waits for a round of
# other clones. After each job, a small checkpoint of that clone is saved for 
# resumability. 
#
# Layout of the checkpoint folder: <dir>/checkpoint_<id>/
#   status.tsv              one row per clone: status and number of rounds done
#   logs/                   output of the worker processes (parallel runs only)
#   clone_<clone_id>/       last good checkpoint of <clone_id>:
#   |_ manifest.tsv         details below 
#   |_result.rds            one-row tibble of the clone
#   |_*.xml, *.state       copies of the small BEAST files 
#
# .log, .trees and the console log are not copied: BEAST just adds to them,
# so we can just keep track of their size at the checkpoint, and cut them to that
# size when we "roll back" to the checkpoint. The small files are copied back. 
#
# Note: tree logs ends with "End;", and BEAST replaces it when it resumes, 
# so we calculate size without the "End;".

# determine which clones to start from beginning due to checkpoint options or
# XML not matching
getClonesToRestart <- function(dir, id, clone_ids, clone_options, xml_md5){
  cp_dir <- getCheckpointDir(dir, id)
  changed <- vapply(clone_ids, function(cid){
    m <- parseIterateManifest(file.path(cp_dir, paste0("clone_", cid), "manifest.tsv"))
    if(is.null(m)){
      return(FALSE)
    }
    old_xml <- m$md5[match(paste0(id, "_", cid, ".xml"), m$files)]
    if(is.na(old_xml)){
      return(FALSE)
    }
    if(!identical(m$options, clone_options(cid))){
      return(TRUE)
    }
    new_xml <- xml_md5(cid)
    !is.na(new_xml) && !identical(new_xml, old_xml)
  }, logical(1))
  restart <- clone_ids[changed]
  if(length(restart) > 0){
    message("The seed, the java flag, BEAST or the XML (the data, the template or an argument) of clone(s) ",
      paste(restart, collapse=", "), " changed since the earlier run. Their checkpoints and ",
      "BEAST files are deleted, and they start again from round 1. The other clones continue.")
  }
  restart
}

# Seconds as "4m10s" ("-" below 0.5).
formatSeconds <- function(x){
  s <- as.integer(round(x))
  ifelse(is.na(x) | x < 0.5, "-", sprintf("%dm%02ds", s %/% 60L, s %% 60L))
}

# Table printed at the end: clone_id, status, rounds, below_ESS, time.
iterateSummaryTable <- function(clone_ids, status, rounds, state, ignore,
    ess_cutoff, elapsed){
  below <- vapply(clone_ids, function(cid){
    # a clone that has not run has no parameters column
    if(rounds[[cid]] == 0 || !"parameters" %in% names(state[[cid]])){
      return(NA_integer_)
    }
    params <- state[[cid]][["parameters"]][[1]]
    tryCatch(as.integer(countParametersBelowESS(params, ignore, ess_cutoff)),
      error=function(e)NA_integer_)
  }, integer(1))
  data.frame(clone_id=clone_ids, status=unname(status[clone_ids]),
    rounds=unname(rounds[clone_ids]), below_ESS=unname(below),
    time=formatSeconds(unname(elapsed[clone_ids])),
    stringsAsFactors=FALSE)
}

# Path of the checkpoint folder.
getCheckpointDir <- function(dir, id){
  file.path(dir, paste0("checkpoint_", id))
}

# Number of parameters with ESS below ess_cutoff, without those that match `ignore`.
countParametersBelowESS <- function(params, ignore, ess_cutoff){
  if(!is.data.frame(params) || !all(c("item", "ESS") %in% names(params))){
    stop("parameter table has no item and ESS columns")
  }
  for(regex in ignore){
    params <- params[!grepl(regex, params$item), , drop=FALSE]
  }
  sum(params$ESS[!params$item %in% ignore] < ess_cutoff, na.rm=TRUE)
}

# Get files in `dir` that belong to each clone.  
getCloneFiles <- function(dir, id, clone_ids){
  files <- list.files(dir)
  files <- files[!dir.exists(file.path(dir, files))] # remove subfolders
  # File names for clone "c": <id>_<c>.<ext> or <id>_<c>_<more>.<ext>
  prefix <- paste0(id, "_", clone_ids)
  owner <- rep(NA_integer_, length(files))
  # Clone ids can be "a" and "a_b" so use longest matching prefix
  for(i in order(nchar(prefix), decreasing=TRUE)){      # longest prefix first
    n <- nchar(prefix[i])
    hit <- is.na(owner) & startsWith(files, prefix[i]) &
      substr(files, n + 1, n + 1) %in% c(".", "_")
    owner[hit] <- i
  }
  out <- split(files, factor(owner, levels=seq_along(clone_ids)))
  names(out) <- clone_ids
  out
}

# Write a text file (write, then rename), so a crash leaves no half-written file.
writeLinesSafely <- function(lines, path){
  tmp <- paste0(path, ".tmp")
  writeLines(lines, tmp)
  if(!file.rename(tmp, path)){
    file.copy(tmp, path, overwrite=TRUE)
    unlink(tmp)
  }
  invisible(path)
}

# Size to save/cut a growing file at the end of a round. 
# BEAST replaces the closing "End;" of tree logs on resume, so save size without it.
getCutSize <- function(path){
  size <- file.size(path)
  n <- min(size, 16)
  con <- file(path, "rb")
  seek(con, size - n)
  tail_bytes <- readBin(con, "raw", n)
  close(con)
  text <- rawToChar(tail_bytes[tail_bytes != as.raw(0)])
  pos <- as.integer(regexpr("End;[ \t\r\n]*$", text, useBytes=TRUE))
  if(pos > 0) size - (nchar(text, type="bytes") - pos + 1) else size
}

# Copy the first `size` bytes of a file to a new file, `chunk` bytes at a time.
iterateCopyPrefix <- function(from, to, size, chunk=1e7){
  inp <- file(from, "rb")
  out <- file(to, "wb")
  left <- size
  while(left > 0){
    bytes <- readBin(inp, "raw", min(left, chunk))
    if(length(bytes) == 0){
      break
    }
    writeBin(bytes, out)
    left <- left - length(bytes)
  }
  close(inp)
  close(out)
  invisible(TRUE)
}

# Cut a file to `size` bytes. 
# Handles the case where truncate() does not work on this system.
# For tests, use use_truncate=FALSE.
truncateFile <- function(path, size, use_truncate=TRUE, chunk=1e7){
  # Nothing to cut if file ends at `size` or size + "End;".
  if(getCutSize(path) == size){
    return(invisible(TRUE))
  }
  if(use_truncate){
    con <- file(path, "r+b")
    ok <- tryCatch({
      seek(con, size, rw="write")
      truncate(con)
      TRUE
    }, error=function(e)FALSE)
    close(con)
    if(ok){
      return(invisible(TRUE))
    }
  }
  tmp <- paste0(path, ".cut")
  iterateCopyPrefix(path, tmp, size, chunk)
  if(!file.rename(tmp, path)){
    file.copy(tmp, path, overwrite=TRUE)
    unlink(tmp)
  }
  invisible(TRUE)
}

# Save the checkpoint of a clone after `round` rounds. 
# `options` are the non-XML run parameters. 
# Should be called only when there are no running BEAST job runs for the clone.
# Small files are copied, size is saved for growing files. 
# Built in clone_<id>.new and swapped in, so there is always an intact checkpoint. 
writeIterateCheckpoint <- function(dir, id, clone_ids, clone_id, round, result, options=""){
  cp_dir <- getCheckpointDir(dir, id)
  files <- getCloneFiles(dir, id, clone_ids)[[clone_id]]
  files <- files[!(grepl("\\.tree$", files) | grepl("_log\\.tsv$", files))]
  grows <- grepl("\\.(log|trees)$", files)
  final <- file.path(cp_dir, paste0("clone_", clone_id))
  tmp <- paste0(final, ".new")
  unlink(tmp, recursive=TRUE)
  if(!dir.create(tmp, recursive=TRUE)){
    stop("Could not create checkpoint folder ", tmp)
  }
  if(!all(file.copy(file.path(dir, files[!grows]), tmp))){
    stop("Could not copy all files of clone ", clone_id, " to ", tmp)
  }
  saveRDS(result, file.path(tmp, "result.rds"))
  copies <- c(files[!grows], "result.rds")
  sizes <- c(vapply(file.path(dir, files[grows]), getCutSize, numeric(1), USE.NAMES=FALSE),
    file.size(file.path(tmp, copies)))
  md5 <- c(rep("-", sum(grows)), unname(tools::md5sum(file.path(tmp, copies))))
  kind <- rep(c("append", "copy"), c(sum(grows), length(copies)))
  writeLinesSafely(c(paste("round", round, sep="\t"), paste0("options\t", options),
    paste(c(files[grows], copies), kind, sprintf("%.0f", sizes), md5, sep="\t")),
    file.path(tmp, "manifest.tsv"))
  unlink(final, recursive=TRUE)
  if(!file.rename(tmp, final)){
    stop("Could not move checkpoint folder ", tmp, " to ", final)
  }
  invisible(final)
}

# Read manifest.tsv: list(round, options, files, kind, size, md5). 
# One value per file. kind is "append" (growing file: only save size) or 
# "copy" (copy in checkpoint folder).
# NULL if the file is missing or not valid.
parseIterateManifest <- function(manifest){
  if(!file.exists(manifest)){
    return(NULL)
  }
  lines <- readLines(manifest, warn=FALSE)
  if(length(lines) < 3 || !grepl("^round\t[0-9]+$", lines[1]) ||
      !startsWith(lines[2], "options\t")){
    return(NULL)
  }
  parts <- strsplit(lines[-(1:2)], "\t", fixed=TRUE)
  if(any(lengths(parts) != 4)){
    return(NULL)
  }
  field <- function(i) vapply(parts, function(x)x[i], character(1))
  size <- suppressWarnings(as.numeric(field(3)))
  if(anyNA(size) || !all(field(2) %in% c("append", "copy"))){
    return(NULL)
  }
  list(round=as.integer(sub("^round\t", "", lines[1])), options=substring(lines[2], 9),
    files=field(1), kind=field(2), size=size, md5=field(4))
}

# Read the last good checkpoint of a clone.
# Returns its manifest (see parseIterateManifest), its checkpoint folder and the
# result of its last round, or NULL if no good checkpoint exists.
# If clone_<id> is not good, clone_<id>.new is tried next.
readCheckpoint <- function(cp_dir, clone_id, dir=NULL){
  final <- file.path(cp_dir, paste0("clone_", clone_id))
  for(cand in c(final, paste0(final, ".new"))){
    m <- parseIterateManifest(file.path(cand, "manifest.tsv"))
    if(is.null(m) || !"result.rds" %in% m$files){
      next
    }
    copies <- m$files[m$kind == "copy"]
    paths <- file.path(cand, copies)
    if(!all(file.exists(paths)) ||
        !isTRUE(all(unname(tools::md5sum(paths)) == m$md5[match(copies, m$files)]))){
      next
    }
    if(!is.null(dir)){
      grown <- m$files[m$kind == "append"]
      paths <- file.path(dir, grown)
      if(!all(file.exists(paths)) || any(file.size(paths) < m$size[match(grown, m$files)])){
        next
      }
    }
    result <- tryCatch(readRDS(file.path(cand, "result.rds")), error=function(e)NULL)
    if(is.null(result)){
      next
    }
    return(c(m, list(dir=cand, result=result)))
  }
  NULL
}

# Put the files of a clone back to the state of its checkpoint `cp`:
# files not in `cp` are deleted, growing files (tree and log) are cut to `cp` size,
# small files are copied back.
doCheckpointRollback <- function(dir, id, clone_ids, clone_id, cp, use_truncate=TRUE){
  owned <- getCloneFiles(dir, id, clone_ids)[[clone_id]]
  owned <- owned[!(grepl("\\.tree$", owned) | grepl("_log\\.tsv$", owned))]
  unlink(file.path(dir, setdiff(owned, cp$files)))
  files <- setdiff(cp$files, "result.rds")
  kind <- cp$kind[match(files, cp$files)]
  for(f in files[kind == "append"]){
    truncateFile(file.path(dir, f), cp$size[match(f, cp$files)], use_truncate)
  }
  if(!all(file.copy(file.path(cp$dir, files[kind == "copy"]), dir, overwrite=TRUE))){
    stop("Could not restore all files of clone ", clone_id, " from ", cp$dir)
  }
  invisible(NULL)
}

# Arguments for one clone: a function in `dots` is called with the clone's row,
# and `trees` is cut to the tree of the clone.
iterateCloneDots <- function(dots, cs, all_ids){
  out <- lapply(dots, function(d) if(is.function(d)) d(cs) else d)
  if(!is.null(dots[["trees"]]) && !is.function(dots[["trees"]])){
    out[["trees"]] <- dots[["trees"]][match(as.character(cs$clone_id), all_ids)]
  }
  out
}

# An error with only its message (its call can hold all the data of a clone).
iterateSlimError <- function(e){
  simpleError(conditionMessage(e))
}

# Runs the jobs one after the other in this R process.
iterateExecutorSync <- function(){
  pending <- list()
  list(
    submit=function(key, work, log_file){
      pending[key] <<- list(tryCatch(work(), error=iterateSlimError))
      invisible(NULL)
    },
    collect=function(){
      out <- pending
      pending <<- list()
      out
    },
    n_active=function() length(pending),
    cleanup=function() invisible(NULL)
  )
}

# Runs each job in a forked R process (not on Windows). Its output goes to its own log file.
iterateExecutorFork <- function(){
  jobs <- list()   # names are process ids
  list(
    submit=function(key, work, log_file){
      job <- parallel::mcparallel({
        con <- file(log_file, open="a")
        sink(con, append=TRUE)
        sink(con, append=TRUE, type="message")
        out <- tryCatch(work(), error=iterateSlimError)
        sink(type="message")
        sink()
        close(con)
        out
      })
      jobs[[as.character(job$pid)]] <<- list(job=job, key=key)
      invisible(NULL)
    },
    collect=function(){
      if(length(jobs) == 0){
        return(list())
      }
      res <- parallel::mccollect(lapply(jobs, function(x)x$job),
        wait=FALSE, timeout=1)
      if(is.null(res)){
        return(list())
      }
      out <- list()
      for(pid in names(res)){
        key <- jobs[[pid]]$key
        value <- res[[pid]]
        if(is.null(value) || inherits(value, "try-error")){
          value <- simpleError("worker process returned no result (it may have been stopped)")
        }
        out[key] <- list(value)
        jobs[[pid]] <<- NULL
      }
      out
    },
    n_active=function() length(jobs),
    cleanup=function(){
      for(j in jobs){
        tryCatch(tools::pskill(j$job$pid, tools::SIGTERM), error=function(e)NULL)
      }
      if(length(jobs) > 0){
        tryCatch(parallel::mccollect(lapply(jobs, function(x)x$job),
          wait=TRUE, timeout=5), error=function(e)NULL)
      }
      jobs <<- list()
      invisible(NULL)
    }
  )
}

# Runs all clones, one BEAST job per clone, until each one converges or has run `iterations`
# rounds, and saves a checkpoint after each round. 
# Returns the clones that did not fail, with below_ESS.
runTimeTreesIterate <- function(clones, work, dir, id, iterations=10, ess_cutoff=200,
    ignore=c("traitfrequencies"), continue=FALSE, nproc=1, quiet=0,
    skip_errors=FALSE, describe=NULL, beast_options=NULL){

  t_start <- Sys.time()

  clone_ids <- as.character(clones$clone_id)
  n <- length(clone_ids)
  if(n == 0){
    clones$below_ESS <- integer(0)
    return(clones)
  }
  if(anyDuplicated(clone_ids) > 0){
    stop("clone_id values must be unique")
  }
  nproc <- max(1, nproc)
  use_fork <- nproc > 1 && .Platform$OS.type != "windows"
  if(nproc > 1 && !use_fork){
    warning("Running clones in parallel needs forked processes, which Windows does not have. Using nproc=1.")
    nproc <- 1
  }

  dir <- path.expand(dir)
  if(!dir.exists(dir)){
    dir.create(dir, recursive=TRUE)
  }
  cp_dir <- getCheckpointDir(dir, id)
  status_path <- file.path(cp_dir, "status.tsv")
  log_dir <- file.path(cp_dir, "logs")

  state <- stats::setNames(lapply(seq_len(n), function(i)clones[i, ]), clone_ids)
  rounds <- stats::setNames(integer(n), clone_ids)
  status <- stats::setNames(rep("queued", n), clone_ids)
  # seconds that each clone used in this call, and the time of its last start
  elapsed <- stats::setNames(numeric(n), clone_ids)
  started <- stats::setNames(rep(NA_real_, n), clone_ids)
  # TRUE after the describe text of a clone was printed
  described <- stats::setNames(rep(FALSE, n), clone_ids)

  if(isTRUE(continue)){
    old_dir <- file.path(dir, paste0("iterations_", id))
    if(!dir.exists(cp_dir) && dir.exists(old_dir)){
      warning("Found ", old_dir, " from an older version of dowser. ",
        "It cannot be used. All clones will start again from round 1.")
    }
    for(cid in clone_ids){
      cp <- readCheckpoint(cp_dir, cid, dir)
      if(is.null(cp)){
        final <- file.path(cp_dir, paste0("clone_", cid))
        if(dir.exists(final) || dir.exists(paste0(final, ".new"))){
          message("The checkpoint of clone ", cid,
            " is damaged. This clone will start again from round 1.")
        }
        next
      }
      doCheckpointRollback(dir, id, clone_ids, cid, cp)
      state[[cid]] <- cp$result
      # Use the data of this call. The stored row has the data of the earlier
      # call. If the data changed, the new XML is not the same as the XML in the
      # checkpoint, and the XML check stops the clone before BEAST starts.
      if("data" %in% names(clones) && "data" %in% names(state[[cid]])){
        state[[cid]][["data"]] <- clones[["data"]][match(cid, clone_ids)]
      }
      rounds[cid] <- cp$round
    }
  }else{
    if(dir.exists(cp_dir)){
      warning("Deleting the existing checkpoints in ", cp_dir,
        ". To continue an earlier run, use continue=TRUE.")
      unlink(cp_dir, recursive=TRUE)
    }
  }
  dir.create(log_dir, recursive=TRUE, showWarnings=FALSE)

  # decide which clones need a round
  queue <- character(0)
  for(cid in clone_ids){
    if(rounds[[cid]] == 0){
      queue <- c(queue, cid)
      next
    }
    below <- countParametersBelowESS(state[[cid]]$parameters[[1]], ignore, ess_cutoff)
    if(below > 0 && rounds[[cid]] < iterations){
      queue <- c(queue, cid)
    }else{
      status[cid] <- if(below > 0) "unconverged" else "done"
    }
  }
  failed <- character(0)

  write_status <- function(){
    lines <- c("clone_id\tstatus\trounds_done",
      paste(clone_ids, status[clone_ids], rounds[clone_ids], sep="\t"))
    tryCatch(writeLinesSafely(lines, status_path), error=function(e)NULL)
  }
  print_summary <- function(){
    tbl <- table(factor(status, levels=c("queued", "running", "done",
      "unconverged", "failed")))
    cat("=== status", format(Sys.time()), "===\n")
    cat(paste(names(tbl), as.integer(tbl), sep="=", collapse=" "), "\n")
    running <- clone_ids[status[clone_ids] == "running"]
    if(length(running) > 0){
      cat("running:", paste0(running, "(round ", rounds[running] + 1L, ")",
        collapse=" "), "\n")
    }
  }

  # With skip_errors=FALSE, the first error does not stop the other jobs.
  # No new clone starts. Clones that are running finish their current run and
  # are saved in checkpoints. Then the function stops with the error messages.
  # No process is stopped, so no BEAST process is left running.
  errors <- character(0)
  handle_failure <- function(cid, msg){
    msg <- paste0("Clone ", cid, " failed: ", msg)
    failed <<- c(failed, cid)
    status[cid] <<- "failed"
    if(skip_errors){
      warning(msg, call.=FALSE)
      return(invisible(NULL))
    }
    errors <<- c(errors, msg)
    if(length(errors) == 1){
      queue <<- character(0)
      text <- paste0(msg, "\nNo new clone will start.")
      if(executor$n_active() > 0){
        text <- paste0(text, " Waiting for ", executor$n_active(),
          " running clone(s) to finish their current run.")
      }
      message(text)
    }
    invisible(NULL)
  }

  executor <- if(use_fork) iterateExecutorFork() else iterateExecutorSync()
  on.exit(executor$cleanup(), add=TRUE)

  launch <- function(cid){
    this_round <- rounds[[cid]] + 1L
    resume_id <- if(rounds[[cid]] == 0) NULL else cid
    cs <- state[[cid]]
    if(quiet < 1){
      # In a sequential run, the worker prints its own output to the console
      extra <- ""
      if(use_fork && !is.null(describe) && !described[[cid]]){
        extra <- tryCatch(paste0("(", describe(cs), ")"), error=function(e)"")
        described[cid] <<- TRUE
      }
      cat("Starting clone", cid, "round", this_round, extra, "\n")
    }
    started[cid] <<- as.numeric(Sys.time())
    status[cid] <<- "running"
    write_status()
    job <- function(){
      res <- work(cs, resume_id)
      opts <- if(is.null(beast_options)) "" else beast_options(cs)
      writeIterateCheckpoint(dir, id, clone_ids, cid, this_round, res, opts)
      res
    }
    executor$submit(cid, job,
      log_file=file.path(log_dir, paste0("clone", cid, "_dowser.log")))
  }
  fill_slots <- function(){
    while(length(queue) > 0 && executor$n_active() < nproc){
      launch(queue[1])
      queue <<- queue[-1]
    }
  }

  write_status()
  fill_slots()
  while(executor$n_active() > 0){
    finished <- executor$collect()
    if(length(finished) == 0){
      next
    }
    for(cid in names(finished)){
      result <- finished[[cid]]
      elapsed[cid] <- elapsed[[cid]] + (as.numeric(Sys.time()) - started[[cid]])
      if(inherits(result, "error")){
        handle_failure(cid, conditionMessage(result))
        next
      }
      if(!inherits(result, "data.frame") || nrow(result) != 1 ||
          !"parameters" %in% names(result)){
        handle_failure(cid, "the worker did not return one row with parameters")
        next
      }
      below <- tryCatch(countParametersBelowESS(result$parameters[[1]], ignore, ess_cutoff),
        error=function(e)e)
      if(inherits(below, "error")){
        handle_failure(cid, conditionMessage(below))
        next
      }
      state[[cid]] <- result
      rounds[cid] <- rounds[[cid]] + 1L
      if(below > 0 && rounds[[cid]] < iterations){
        # after an error, the clone is not queued, but continue=TRUE runs it
        if(length(errors) == 0){
          queue <- c(queue, cid)
        }
        status[cid] <- "queued"
      }else{
        status[cid] <- if(below > 0) "unconverged" else "done"
      }
    }
    fill_slots()
    write_status()
    if(quiet < 1){
      print_summary()
    }
  }

  if(quiet < 1){
    cat("=== summary", format(Sys.time()), "===\n")
    print(iterateSummaryTable(clone_ids, status, rounds, state, ignore,
      ess_cutoff, elapsed), row.names=FALSE)
    cat("time in this call:", formatSeconds(
      as.numeric(difftime(Sys.time(), t_start, units="secs"))), "\n")
  }

  if(length(errors) > 0){
    stop(paste(errors, collapse="\n"),
      "\nClones that finished a run are saved in ", cp_dir,
      ". To resume the analysis, run it again with continue=TRUE.", call.=FALSE)
  }

  unconverged <- clone_ids[status[clone_ids] == "unconverged"]
  if(length(unconverged) > 0){
    warning(paste(paste(unconverged, collapse=","),
      "failed to converge after", iterations, "iterations"))
  }
  if(length(failed) > 0){
    warning("Dropped clones due to errors: ", paste(failed, collapse=","))
  }

  keep <- setdiff(clone_ids, failed)
  out <- dplyr::bind_rows(state[keep])
  out$below_ESS <- vapply(out$parameters, countParametersBelowESS, integer(1),
    ignore=ignore, ess_cutoff=ess_cutoff)
  out
}

# Move (or copy if keep=TRUE) a file or folder, or copy it if keep=TRUE. 
moveFiles <- function(from, to, keep=FALSE){
  if(!keep && file.rename(from, to)){
    return(invisible(TRUE))
  }
  if(dir.exists(from)){
    ok <- file.copy(from, dirname(to), recursive=TRUE)
  }else{
    ok <- file.copy(from, to)
  }
  if(!all(ok)){
    stop("Could not ", if(keep) "copy " else "move ", from, " to ", to)
  }
  if(!keep){
    unlink(from, recursive=TRUE)
  }
  invisible(TRUE)
}

# Make backup folder. If continue=TRUE, the files are copied, otherwise they are moved.
makeBackup <- function(dir, id, clone_ids, continue){
  files <- unlist(getCloneFiles(dir, id, clone_ids), use.names=FALSE)
  cp_dir <- getCheckpointDir(dir, id)
  has_cp <- dir.exists(cp_dir)
  if(length(files) == 0 && !has_cp){
    return(NULL)
  }
  backup_name <- paste0("backup_", id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  backup_dir <- file.path(dir, backup_name)
  k <- 1
  while(dir.exists(backup_dir)){
    k <- k + 1
    backup_dir <- file.path(dir, paste0(backup_name, "_", k))
  }
  dir.create(backup_dir, recursive=TRUE)
  for(f in files){
    moveFiles(file.path(dir, f), file.path(backup_dir, f), keep=continue)
  }
  if(has_cp){
    moveFiles(cp_dir, file.path(backup_dir, basename(cp_dir)), keep=continue)
  }
  writeLines(c("Backup made by getTimeTreesIterate(backup=TRUE).",
    paste("id:", id), paste("time:", format(Sys.time())),
    "To put these files back, stop the new run and use restoreIterateBackup()."),
    file.path(backup_dir, "README.txt"))

  if(!is.null(backup_dir)){
    message("The existing files of analysis '", id, "' were ",
      if(isTRUE(continue)) "copied" else "moved", " to ", backup_dir,
      ".\nTo undo this run, stop it and use restoreIterateBackup(dir, id). ",
      "Delete the backup folder yourself when you do not need it any more.")
  }
  backup_dir
}

#' Put back the files that getTimeTreesIterate moved to a backup folder
#'
#' Assumes no run or BEAST processes are running with this \code{dir} and \code{id}.
#' Deletes any new checkpoint files since the last backup. Restores checkpoint and 
#' BEAST files from the last or specified backup. The backup folder is deleted.
#'
#' @param    dir     directory of the analysis
#' @param    id      analysis id
#' @param    backup  backup folder. Default: newest backup of \code{id} in \code{dir}.
#'
#' @return   Invisibly, the names of the files and folders that were put back.
#' @seealso \link{getTimeTreesIterate}
#' @export
restoreIterateBackup <- function(dir, id, backup=NULL){
  dir <- path.expand(dir)
  if(is.null(backup)){
    prefix <- paste0("backup_", id, "_")
    found <- list.files(dir, full.names=FALSE)
    found <- found[startsWith(found, prefix) & dir.exists(file.path(dir, found))]
    found <- found[grepl("^[0-9]{8}_[0-9]{6}(_[0-9]+)?$", substring(found, nchar(prefix) + 1))]
    if(length(found) == 0){
      stop("No backup folder found for id '", id, "' in ", dir)
    }
    backup <- file.path(dir, sort(found)[length(found)])
  }
  backup <- path.expand(backup)
  if(!dir.exists(backup)){
    stop("Backup folder ", backup, " not found")
  }
  items <- setdiff(list.files(backup, all.files=FALSE), "README.txt")
  # the checkpoints of the new run do not fit the old files
  unlink(getCheckpointDir(dir, id), recursive=TRUE)
  for(item in items){
    target <- file.path(dir, item)
    unlink(target, recursive=TRUE)
    moveFiles(file.path(backup, item), target)
  }
  unlink(backup, recursive=TRUE)
  message("Put back ", length(items), " files and folders from ", backup)
  invisible(items)
}

#' Read the results of a run as they were at its last checkpoint
#'
#' Use if a \code{getTimeTreesIterate} run stopped and you want
#' the results as they were at the end of each clone's last finished round. 
#'
#' @param    clones      tibble of \code{airrClone} objects, as in \code{getTimeTreesIterate}
#' @param    dir         directory of the analysis
#' @param    id          analysis id
#' @param    beast       location of the beast binary directory (beast/bin)
#' @param    burnin      percent of initial samples to discard (default 10)
#' @param    trait       trait column that was used, if any
#' @param    nproc       number of cores for parallelization. Uses at most 1 core per clone.
#' @param    quiet       amount of rubbish to print to console
#' @param    posterior   what to read of the full posterior: "none" (default), "all", or a vector with "parameters", "trees_with_traits" and "trees". See \code{readBEAST}.
#' @param    asr         were ancestral sequences logged?
#' @param    low_ram     run with less memory (slightly slower)
#' @param    trim_ids    remove last _ group from tips?
#' @param    ess_cutoff  minimum ESS, for the column below_ESS
#' @param    ignore      vector of parameters to ignore for the column below_ESS
#'
#' @return   Tibble of clones that have a good checkpoint, with additional columns
#' \code{trees}, \code{parameters}, \code{below_ESS} and \code{checkpoint_round} (number of rounds finished for that clone). A warning lists the clones that have no good checkpoint.
#' @seealso \link{getTimeTreesIterate}, \link{readBEAST}
#' @export
readFromLastCheckpoint <- function(clones, dir, id, beast, burnin=10, trait=NULL, nproc=1, quiet=0,
    posterior=c("none", "all", "parameters", "trees_with_traits", "trees"), asr=FALSE,
    low_ram=TRUE, trim_ids=FALSE, ess_cutoff=200, ignore=c("traitfrequencies")){
  dir <- path.expand(dir)
  # the default vector of readBEAST means "all", but the default here is "none"
  if(setequal(posterior, c("none", "all", "parameters", "trees_with_traits", "trees"))){
    posterior <- "none"
  }
  clone_ids <- as.character(clones$clone_id)
  cp_dir <- getCheckpointDir(dir, id)
  cps <- lapply(stats::setNames(clone_ids, clone_ids),
    function(cid) readCheckpoint(cp_dir, cid, dir))
  bad <- clone_ids[vapply(cps, is.null, logical(1))]
  if(length(bad) > 0){
    warning("No good checkpoint for clone(s) ", paste(bad, collapse=", "),
      ". They are not in the result.")
  }
  keep <- setdiff(clone_ids, bad)
  if(length(keep) == 0){
    stop("No clone has a good checkpoint in ", cp_dir)
  }

  # rebuild the files that readBEAST needs, as they were at the end of the last round
  tmp <- tempfile("checkpoint_read_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
  for(cid in keep){
    cp <- cps[[cid]]
    for(f in setdiff(cp$files, "result.rds")){
      i <- match(f, cp$files)
      if(cp$kind[i] == "copy"){
        if(!file.copy(file.path(cp$dir, f), file.path(tmp, f))){
          stop("Could not copy ", f, " from the checkpoint of clone ", cid)
        }
      }else{
        iterateCopyPrefix(file.path(dir, f), file.path(tmp, f), cp$size[i])
        if(grepl("\\.trees$", f)){
          # the closing line that BEAST had written, and removes when it resumes
          cat("End;", file=file.path(tmp, f), append=TRUE)
        }
      }
    }
    if(quiet < 1){
      cat("Clone", cid, ": results as of round", cp$round, "\n")
    }
  }

  out <- readBEAST(clones=clones[match(keep, clone_ids), ], dir=tmp, id=id, beast=beast,
    burnin=burnin, trait=trait, nproc=nproc, quiet=quiet, posterior=posterior, asr=asr,
    low_ram=low_ram, trim_ids=trim_ids)
  out$checkpoint_round <- vapply(keep, function(cid) cps[[cid]]$round, integer(1), USE.NAMES=FALSE)
  out$below_ESS <- vapply(out$parameters, countParametersBelowESS, integer(1),
    ignore=ignore, ess_cutoff=ess_cutoff)
  out
}
