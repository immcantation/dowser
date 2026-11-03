# Tests for getTimeTreesIterate and what goes with it (R/TimeTreesIterateCore.R, and the wrapper in
# R/TimeTreesFunctions.R). BEAST is never run.
#
# Part 1: the helper functions. In each round the fake worker adds one line to three fake BEAST
#         files of the clone (.log, .trees, .xml.state). The number of lines in the .log file is
#         the number of the round.
# Part 2: the engine (queue, checkpoints, rollback, truncate) and what continue=TRUE can use.
# Part 3: the wrapper getTimeTreesIterate, restoreIterateBackup and readFromLastCheckpoint, with
#         the backup helpers. getTimeTrees and readBEAST are replaced by fake functions. These
#         tests need testthat 3.2.0 or later (local_mocked_bindings).

# ============================================================================================
# Part 1: helper functions
# ============================================================================================

fake_params <- function(ess){
  data.frame(item=c("posterior", "likelihood", "traitfrequencies.1"),
    ESS=c(ess, ess, 1), stringsAsFactors=FALSE)
}

# ess_fn(clone_id, round) gives the ESS. NA makes the worker fail.
iterate_setup <- function(ids, ess_fn){
  dir <- tempfile("iterate_test_")
  dir.create(dir)
  calls <- new.env()
  calls$log <- list()
  work <- function(cs, resume_id){
    cid <- as.character(cs$clone_id)
    for(ext in c(".log", ".trees", ".xml.state")){
      cat("round\n", file=file.path(dir, paste0("t1_", cid, ext)), append=TRUE)
    }
    n <- length(readLines(file.path(dir, paste0("t1_", cid, ".log"))))
    calls$log[[length(calls$log) + 1]] <- list(clone=cid, resume=resume_id, n=n,
      data=cs$data[[1]])
    ess <- ess_fn(cid, n)
    if(is.na(ess)){
      stop("fake failure")
    }
    cs$parameters <- list(fake_params(ess))
    cs
  }
  list(dir=dir, calls=calls, work=work,
    clones=dplyr::tibble(clone_id=ids, data=as.list(ids)))
}

iterate_run <- function(s, ...){
  args <- list(...)
  if(is.null(args$quiet)) args$quiet <- 1
  do.call(dowser:::runTimeTreesIterate, c(list(s$clones, s$work, dir=s$dir, id="t1"), args))
}

# run f and return the value and all warning messages
with_warnings <- function(f){
  msgs <- character(0)
  value <- withCallingHandlers(f(), warning=function(w){
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value=value, warnings=msgs)
}

calls_of <- function(s) s$calls$log
call_clones <- function(s) vapply(s$calls$log, function(x)x$clone, character(1))

# ============================================================================================
# Part 2: the engine, and what continue=TRUE can use
# ============================================================================================

test_that("countParametersBelowESS counts parameters below the cutoff and ignores patterns", {
  expect_equal(dowser:::countParametersBelowESS(fake_params(100), "traitfrequencies", 200), 2)
  expect_equal(dowser:::countParametersBelowESS(fake_params(200), "traitfrequencies", 200), 0)
  # the ignored parameter has ESS 1, so it is counted if it is not ignored
  expect_equal(dowser:::countParametersBelowESS(fake_params(200), character(0), 200), 1)
  expect_error(dowser:::countParametersBelowESS(data.frame(x=1), "a", 200), "ESS")
})

test_that("engine runs until convergence and adds below_ESS", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  out <- iterate_run(s)
  expect_equal(out$clone_id, c("a", "b"))
  expect_equal(out$below_ESS, c(0L, 0L))
  expect_equal(sum(call_clones(s) == "a"), 2)
  expect_equal(sum(call_clones(s) == "b"), 2)
  # first round starts a new chain, later rounds resume it
  for(x in calls_of(s)){
    if(x$n == 1) expect_null(x$resume) else expect_equal(x$resume, x$clone)
  }
})

test_that("engine stops at the iteration limit and warns", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 1)
  res <- with_warnings(function() iterate_run(s, iterations=3))
  expect_true(any(grepl("failed to converge after 3 iterations", res$warnings)))
  expect_equal(res$value$below_ESS, c(2L, 2L))
  expect_equal(length(calls_of(s)), 6)
})

test_that("engine runs a clone only once if it converged in round 1", {
  s <- iterate_setup("a", function(cid, n) 1000)
  out <- iterate_run(s)
  expect_equal(out$below_ESS, 0L)
  expect_equal(length(calls_of(s)), 1)
})

test_that("engine returns the clones with below_ESS when there are no clones", {
  s <- iterate_setup(character(0), function(cid, n) 1000)
  out <- iterate_run(s)
  expect_equal(nrow(out), 0)
  expect_true("below_ESS" %in% names(out))
})

test_that("engine refuses duplicate clone ids", {
  s <- iterate_setup(c("a", "a"), function(cid, n) 1000)
  expect_error(iterate_run(s), "unique")
})

test_that("engine writes checkpoints and a status file", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s)
  cp <- file.path(s$dir, "checkpoint_t1")
  expect_true(file.exists(file.path(cp, "clone_a", "manifest.tsv")))
  expect_true(file.exists(file.path(cp, "clone_a", "result.rds")))
  status <- readLines(file.path(cp, "status.tsv"))
  expect_equal(status, c("clone_id\tstatus\trounds_done", "a\tdone\t2", "b\tdone\t2"))
})

test_that("a checkpoint has copies of the small files, and only the sizes of the growing files", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s)
  clone_a <- file.path(s$dir, "checkpoint_t1", "clone_a")
  expect_equal(sort(list.files(clone_a)), c("manifest.tsv", "result.rds", "t1_a.xml.state"))
  m <- readLines(file.path(clone_a, "manifest.tsv"))
  expect_equal(m[1], "round\t2")
  # two lines of "round" are 12 bytes
  expect_true("t1_a.log\tappend\t12\t-" %in% m)
  expect_true("t1_a.trees\tappend\t12\t-" %in% m)
  expect_true(any(grepl("^t1_a.xml.state\tcopy\t12\t[0-9a-f]{32}$", m)))
  # the files in dir are not touched
  expect_equal(length(readLines(file.path(s$dir, "t1_a.log"))), 2)
  expect_true(file.exists(file.path(s$dir, "t1_a.xml.state")))
})

test_that("a run without continue deletes old checkpoints and warns", {
  s <- iterate_setup("a", function(cid, n) 1000)
  iterate_run(s)
  res <- with_warnings(function() iterate_run(s))
  expect_true(any(grepl("Deleting the existing checkpoints", res$warnings)))
})

test_that("continue restores the last checkpoint before the next round", {
  fail <- FALSE
  ess_fn <- function(cid, n) if(cid == "a" && n == 2 && fail) NA else 1
  s <- iterate_setup(c("a", "b"), ess_fn)
  suppressWarnings(iterate_run(s, iterations=1))
  # Round 2 of clone a fails after the worker changed its files. The files in
  # dir are then not the same as the checkpoint.
  fail <- TRUE
  expect_error(suppressMessages(iterate_run(s, iterations=2, continue=TRUE)),
    "Clone a failed")
  expect_equal(length(readLines(file.path(s$dir, "t1_a.log"))), 2)
  fail <- FALSE
  s$calls$log <- list()
  suppressWarnings(iterate_run(s, iterations=2, continue=TRUE))
  a_call <- Filter(function(x)x$clone == "a", calls_of(s))[[1]]
  # the extra line is gone, so this is round 2 and not round 3
  expect_equal(a_call$n, 2)
  expect_equal(a_call$resume, "a")
  expect_equal(readLines(file.path(s$dir, "t1_a.log")), c("round", "round"))
  expect_true(file.exists(file.path(s$dir, "t1_a.xml.state")))
})

test_that("continue replaces files that do not match each other", {
  fail <- FALSE
  ess_fn <- function(cid, n) if(n == 2 && fail) NA else 1
  s <- iterate_setup("a", ess_fn)
  suppressWarnings(iterate_run(s, iterations=1))
  fail <- TRUE
  expect_error(suppressMessages(iterate_run(s, iterations=2, continue=TRUE)),
    "Clone a failed")
  # the state file is from round 1, but the trees file is from round 2
  writeLines("round", file.path(s$dir, "t1_a.xml.state"))
  cat("round\n", file=file.path(s$dir, "t1_a.trees"), append=TRUE)
  fail <- FALSE
  suppressWarnings(iterate_run(s, iterations=2, continue=TRUE))
  expect_equal(readLines(file.path(s$dir, "t1_a.trees")), c("round", "round"))
  expect_equal(readLines(file.path(s$dir, "t1_a.xml.state")), c("round", "round"))
})

test_that("continue does not change the files of clones that ended well", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s)
  before <- tools::md5sum(file.path(s$dir, c("t1_a.log", "t1_b.trees")))
  s$calls$log <- list()
  out <- iterate_run(s, continue=TRUE)
  expect_equal(length(calls_of(s)), 0)
  expect_equal(out$below_ESS, c(0L, 0L))
  expect_equal(tools::md5sum(file.path(s$dir, c("t1_a.log", "t1_b.trees"))), before)
})

test_that("continue starts a clone again if a growing file is shorter than at its checkpoint", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s)
  writeLines("x", file.path(s$dir, "t1_a.log"))   # 2 bytes, and it was 12 at the checkpoint
  s$calls$log <- list()
  expect_message(iterate_run(s, continue=TRUE), "checkpoint of clone a is damaged")
  expect_equal(call_clones(s)[1], "a")
  expect_null(calls_of(s)[[1]]$resume)
})

test_that("continue starts a clone again if a growing file is missing", {
  s <- iterate_setup("a", function(cid, n) 100 * n)
  iterate_run(s)
  unlink(file.path(s$dir, "t1_a.trees"))
  s$calls$log <- list()
  expect_message(iterate_run(s, continue=TRUE), "is damaged")
  expect_null(calls_of(s)[[1]]$resume)
})

test_that("continue cuts the growing files back and deletes the files of a round that did not end", {
  fail <- FALSE
  ess_fn <- function(cid, n) if(n == 2 && fail) NA else 1
  s <- iterate_setup("a", ess_fn)
  suppressWarnings(iterate_run(s, iterations=1))
  fail <- TRUE
  expect_error(suppressMessages(iterate_run(s, iterations=2, continue=TRUE)), "Clone a failed")
  expect_equal(length(readLines(file.path(s$dir, "t1_a.log"))), 2)
  # a file that the failed round made
  writeLines("partial", file.path(s$dir, "t1_a_extra.log"))
  fail <- FALSE
  s$calls$log <- list()
  suppressWarnings(iterate_run(s, iterations=2, continue=TRUE))
  expect_false(file.exists(file.path(s$dir, "t1_a_extra.log")))
  expect_equal(readLines(file.path(s$dir, "t1_a.log")), c("round", "round"))
})

test_that("continue with a higher iterations limit runs a finished clone again", {
  s <- iterate_setup("a", function(cid, n) if(n >= 3) 1000 else 1)
  suppressWarnings(iterate_run(s, iterations=1))
  s$calls$log <- list()
  out <- iterate_run(s, iterations=5, continue=TRUE)
  expect_equal(vapply(calls_of(s), function(x)x$n, numeric(1)), c(2, 3))
  expect_equal(out$below_ESS, 0L)
})

test_that("continue does not run clones that already converged", {
  ess_fn <- function(cid, n) if(cid == "a") 1000 else 1
  s <- iterate_setup(c("a", "b"), ess_fn)
  suppressWarnings(iterate_run(s, iterations=1))
  s$calls$log <- list()
  out <- suppressWarnings(iterate_run(s, iterations=3, continue=TRUE))
  expect_equal(call_clones(s), c("b", "b"))
  expect_equal(out$clone_id, c("a", "b"))
  expect_equal(out$below_ESS, c(0L, 2L))
})

test_that("continue starts a clone again if its checkpoint is damaged", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 1)
  suppressWarnings(iterate_run(s, iterations=1))
  cat("damage\n", file=file.path(s$dir, "checkpoint_t1", "clone_a", "result.rds"),
    append=TRUE)
  s$calls$log <- list()
  expect_message(suppressWarnings(iterate_run(s, iterations=2, continue=TRUE)),
    "checkpoint of clone a is damaged")
  a_call <- Filter(function(x)x$clone == "a", calls_of(s))[[1]]
  expect_null(a_call$resume)
  b_call <- Filter(function(x)x$clone == "b", calls_of(s))[[1]]
  expect_equal(b_call$resume, "b")
})

test_that("continue without checkpoints runs all clones from the start", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 1000)
  out <- iterate_run(s, continue=TRUE)
  expect_equal(out$below_ESS, c(0L, 0L))
  expect_true(all(vapply(calls_of(s), function(x)is.null(x$resume), logical(1))))
})

test_that("a failed clone stops the run unless skip_errors is TRUE", {
  ess_fn <- function(cid, n) if(cid == "b") NA else 1000
  s <- iterate_setup(c("a", "b"), ess_fn)
  expect_error(suppressMessages(iterate_run(s)), "Clone b failed: fake failure")

  s <- iterate_setup(c("a", "b"), ess_fn)
  res <- with_warnings(function() iterate_run(s, skip_errors=TRUE))
  expect_true(any(grepl("Clone b failed", res$warnings)))
  expect_true(any(grepl("Dropped clones due to errors: b", res$warnings)))
  expect_equal(res$value$clone_id, "a")
  expect_equal(res$value$below_ESS, 0L)
  status <- readLines(file.path(s$dir, "checkpoint_t1", "status.tsv"))
  expect_true("b\tfailed\t0" %in% status)
})

test_that("a clone that failed can be run again with continue", {
  fail <- TRUE
  ess_fn <- function(cid, n) if(cid == "b" && fail) NA else 1000
  s <- iterate_setup(c("a", "b"), ess_fn)
  suppressWarnings(iterate_run(s, skip_errors=TRUE))
  fail <- FALSE
  s$calls$log <- list()
  out <- iterate_run(s, continue=TRUE)
  expect_equal(call_clones(s), "b")
  expect_equal(out$clone_id, c("a", "b"))
})

test_that("parallel run gives the same result as a sequential run", {
  skip_on_os("windows")
  skip_on_cran()
  s <- iterate_setup(c("a", "b", "c"), function(cid, n) 100 * n)
  out <- iterate_run(s, nproc=2)
  expect_equal(out$clone_id, c("a", "b", "c"))
  expect_equal(out$below_ESS, c(0L, 0L, 0L))
  for(cid in c("a", "b", "c")){
    expect_equal(length(readLines(file.path(s$dir, paste0("t1_", cid, ".log")))), 2)
  }
})

test_that("parallel run reports a failed clone", {
  skip_on_os("windows")
  skip_on_cran()
  s <- iterate_setup(c("a", "b"), function(cid, n) if(cid == "b") NA else 1000)
  expect_error(suppressMessages(iterate_run(s, nproc=2)), "Clone b failed: fake failure")
})

test_that("getCloneFiles gives each file to the clone with the longest matching name", {
  dir <- tempfile("assign_test_")
  dir.create(dir)
  for(f in c("x_1.log", "x_1.xml.state", "x_1_2.log", "x_1_2.trees", "x_3.log",
      "x_10.log", "y_1.log", "other.txt")){
    file.create(file.path(dir, f))
  }
  dir.create(file.path(dir, "x_1_folder"))
  out <- dowser:::getCloneFiles(dir, "x", c("1", "1_2", "3"))
  expect_equal(sort(out[["1"]]), c("x_1.log", "x_1.xml.state"))
  expect_equal(sort(out[["1_2"]]), c("x_1_2.log", "x_1_2.trees"))
  expect_equal(out[["3"]], "x_3.log")
})

test_that("getCutSize is the size of a file, or the position of a closing End;", {
  write_text <- function(text){
    path <- tempfile("cut_test_")
    writeBin(charToRaw(text), path)
    path
  }
  size_of <- function(text) dowser:::getCutSize(write_text(text))
  expect_equal(size_of("line\nline\n"), 10)
  expect_equal(size_of(""), 0)
  # a tree log: BEAST replaces the closing End; when it resumes
  expect_equal(size_of("tree 1;\nEnd;"), 8)
  expect_equal(size_of("tree 1;\nEnd;\n"), 8)
  expect_equal(size_of("End;"), 0)
  long <- paste0(strrep("0123456789", 100), "):0.0;\nEnd;")
  expect_equal(size_of(long), nchar(long) - 4)
  # End; that is not at the end of the file is not a closing line
  text <- "End; and then more text in the file\n"
  expect_equal(size_of(text), nchar(text))
  text <- paste0("End;", strrep("x", 30))
  expect_equal(size_of(text), nchar(text))
})

test_that("truncateFile cuts a file, with truncate() and without it", {
  for(use in c(TRUE, FALSE)){
    path <- tempfile("truncate_test_")
    writeBin(charToRaw("0123456789"), path)
    dowser:::truncateFile(path, 4, use_truncate=use)
    expect_equal(rawToChar(readBin(path, "raw", 100)), "0123", info=use)
    expect_equal(file.size(path), 4, info=use)
    # the same size: nothing changes
    dowser:::truncateFile(path, 4, use_truncate=use)
    expect_equal(file.size(path), 4, info=use)
    expect_false(file.exists(paste0(path, ".cut")), info=use)
    # to zero bytes
    dowser:::truncateFile(path, 0, use_truncate=use)
    expect_equal(file.size(path), 0, info=use)
  }
})

test_that("truncateFile leaves a tree log that ends with End; after the saved size", {
  for(use in c(TRUE, FALSE)){
    path <- tempfile("truncate_test_")
    writeBin(charToRaw("tree 1;\nEnd;"), path)
    dowser:::truncateFile(path, 8, use_truncate=use)
    expect_equal(rawToChar(readBin(path, "raw", 100)), "tree 1;\nEnd;", info=use)
    # more than End; after the saved size is cut
    writeBin(charToRaw("tree 1;\ntree 2;\nEnd;"), path)
    dowser:::truncateFile(path, 8, use_truncate=use)
    expect_equal(rawToChar(readBin(path, "raw", 100)), "tree 1;\n", info=use)
    # the same for a tree log from a round that did not end (no End;)
    writeBin(charToRaw("tree 1;\ntree 2;\n"), path)
    dowser:::truncateFile(path, 8, use_truncate=use)
    expect_equal(rawToChar(readBin(path, "raw", 100)), "tree 1;\n", info=use)
  }
})

test_that("truncateFile copies a file in more than one chunk", {
  path <- tempfile("truncate_test_")
  writeBin(charToRaw("0123456789"), path)
  dowser:::truncateFile(path, 7, use_truncate=FALSE, chunk=3)
  expect_equal(rawToChar(readBin(path, "raw", 100)), "0123456")
})

test_that("a checkpoint is written and read back", {
  s <- iterate_setup("a", function(cid, n) 1000)
  s$work(s$clones[1, ], NULL)
  result <- s$work(s$clones[1, ], "a")
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 2L, result, "seed=1")
  cp_dir <- dowser:::getCheckpointDir(s$dir, "t1")
  cp <- dowser:::readCheckpoint(cp_dir, "a", s$dir)
  expect_equal(cp$round, 2L)
  expect_equal(cp$options, "seed=1")
  expect_equal(sort(cp$files), c("result.rds", "t1_a.log", "t1_a.trees", "t1_a.xml.state"))
  pos <- function(f) match(f, cp$files)
  expect_equal(cp$kind[pos("t1_a.log")], "append")
  expect_equal(cp$size[pos("t1_a.log")], 12)
  expect_equal(cp$kind[pos("t1_a.xml.state")], "copy")
  expect_equal(cp$kind[pos("result.rds")], "copy")
  expect_equal(cp$result$clone_id, "a")
})

test_that("the size of a tree log in a checkpoint is the position of its End;", {
  s <- iterate_setup("a", function(cid, n) 1000)
  s$work(s$clones[1, ], NULL)
  result <- s$work(s$clones[1, ], "a")
  cat("tree 1;\nEnd;", file=file.path(s$dir, "t1_a.trees"))
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 1L, result)
  cp <- dowser:::readCheckpoint(dowser:::getCheckpointDir(s$dir, "t1"), "a", s$dir)
  expect_equal(cp$size[match("t1_a.trees", cp$files)], 8)
})

test_that("a checkpoint is not used if a growing file is shorter or missing (a longer one is good)", {
  s <- iterate_setup("a", function(cid, n) 1000)
  result <- s$work(s$clones[1, ], NULL)
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 1L, result)
  cp_dir <- dowser:::getCheckpointDir(s$dir, "t1")
  expect_false(is.null(dowser:::readCheckpoint(cp_dir, "a", s$dir)))
  # a growing file that is longer is good: the next round added to it
  cat("more\n", file=file.path(s$dir, "t1_a.log"), append=TRUE)
  expect_false(is.null(dowser:::readCheckpoint(cp_dir, "a", s$dir)))
  # a shorter or a missing growing file is not good
  writeLines("", file.path(s$dir, "t1_a.log"))
  expect_null(dowser:::readCheckpoint(cp_dir, "a", s$dir))
  # ... unless dir is not given, then the growing files are not looked at
  expect_false(is.null(dowser:::readCheckpoint(cp_dir, "a")))
  cat("round\n", file=file.path(s$dir, "t1_a.log"))
  expect_false(is.null(dowser:::readCheckpoint(cp_dir, "a", s$dir)))
  unlink(file.path(s$dir, "t1_a.trees"))
  expect_null(dowser:::readCheckpoint(cp_dir, "a", s$dir))
})

test_that("a checkpoint with a changed copy is not used", {
  s <- iterate_setup("a", function(cid, n) 1000)
  result <- s$work(s$clones[1, ], NULL)
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 1L, result)
  cp_dir <- dowser:::getCheckpointDir(s$dir, "t1")
  cat("x\n", file=file.path(cp_dir, "clone_a", "t1_a.xml.state"), append=TRUE)
  expect_null(dowser:::readCheckpoint(cp_dir, "a"))
})

test_that("a checkpoint with a changed result.rds is not used", {
  s <- iterate_setup("a", function(cid, n) 1000)
  result <- s$work(s$clones[1, ], NULL)
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 1L, result)
  cp_dir <- dowser:::getCheckpointDir(s$dir, "t1")
  cat("x\n", file=file.path(cp_dir, "clone_a", "result.rds"), append=TRUE)
  expect_null(dowser:::readCheckpoint(cp_dir, "a", s$dir))
})

test_that("doCheckpointRollback puts the files of a clone back to the state of its checkpoint", {
  for(use in c(TRUE, FALSE)){
    s <- iterate_setup(c("a", "b"), function(cid, n) 1000)
    s$work(s$clones[1, ], NULL)
    s$work(s$clones[2, ], NULL)
    result <- s$work(s$clones[1, ], "a")      # two lines in each file of clone a
    dowser:::writeIterateCheckpoint(s$dir, "t1", c("a", "b"), "a", 2L, result)
    cp_dir <- dowser:::getCheckpointDir(s$dir, "t1")
    cp <- dowser:::readCheckpoint(cp_dir, "a", s$dir)
    # the next round of clone a adds to the files, and makes a new one. Clone b is not touched.
    s$work(s$clones[1, ], "a")
    writeLines("new", file.path(s$dir, "t1_a_new.log"))
    b_before <- tools::md5sum(file.path(s$dir, "t1_b.log"))
    dowser:::doCheckpointRollback(s$dir, "t1", c("a", "b"), "a", cp, use_truncate=use)
    expect_equal(readLines(file.path(s$dir, "t1_a.log")), c("round", "round"), info=use)
    expect_equal(readLines(file.path(s$dir, "t1_a.trees")), c("round", "round"), info=use)
    expect_equal(readLines(file.path(s$dir, "t1_a.xml.state")), c("round", "round"), info=use)
    expect_false(file.exists(file.path(s$dir, "t1_a_new.log")), info=use)
    expect_equal(tools::md5sum(file.path(s$dir, "t1_b.log")), b_before, info=use)
  }
})

test_that("parseIterateManifest reads a manifest and rejects one that is not correct", {
  path <- tempfile("manifest_")
  good <- c("round\t2", "options\tseed=1 java=TRUE", "a.log\tappend\t12\t-",
    "a.xml\tcopy\t30\tabc", "result.rds\tcopy\t99\tdef")
  writeLines(good, path)
  parsed <- dowser:::parseIterateManifest(path)
  expect_equal(parsed$round, 2L)
  expect_equal(parsed$options, "seed=1 java=TRUE")
  expect_equal(parsed$files, c("a.log", "a.xml", "result.rds"))
  expect_equal(parsed$kind, c("append", "copy", "copy"))
  expect_equal(parsed$size, c(12, 30, 99))
  expect_equal(parsed$md5, c("-", "abc", "def"))
  # an empty options text is valid
  writeLines(c("round\t2", "options\t", "result.rds\tcopy\t9\tabc"), path)
  expect_equal(dowser:::parseIterateManifest(path)$options, "")
  bad <- list(
    c("round\t2", "options\tx"),                                          # no file lines
    c("round\tx", "options\tx", "a.log\tappend\t1\t-"),                  # round is not a number
    c("round\t2", "state\tfull", "options\tx", "a.log\tappend\t1\t-"),   # the format of an older version
    c("round\t2", "options\tx", "a.log\tother\t1\t-"),                   # a kind that does not exist
    c("round\t2", "options\tx", "a.log\tappend\tabc\t-"),                # the size is not a number
    c("round\t2", "options\tx", "a.log\tappend\t1"))                     # a file without md5 field
  for(lines in bad){
    writeLines(lines, path)
    expect_null(dowser:::parseIterateManifest(path), info=paste(lines, collapse=" | "))
  }
  expect_null(dowser:::parseIterateManifest(tempfile("no_manifest_")))
})

test_that("a complete .new checkpoint is used if the old folder is missing", {
  s <- iterate_setup("a", function(cid, n) 1000)
  result <- s$work(s$clones[1, ], NULL)
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 1L, result)
  cp_dir <- dowser:::getCheckpointDir(s$dir, "t1")
  file.rename(file.path(cp_dir, "clone_a"), file.path(cp_dir, "clone_a.new"))
  cp <- dowser:::readCheckpoint(cp_dir, "a")
  expect_false(is.null(cp))
  expect_equal(cp$round, 1L)
})

test_that("a .new checkpoint without manifest is not used", {
  s <- iterate_setup("a", function(cid, n) 1000)
  result <- s$work(s$clones[1, ], NULL)
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 1L, result)
  cp_dir <- dowser:::getCheckpointDir(s$dir, "t1")
  file.rename(file.path(cp_dir, "clone_a"), file.path(cp_dir, "clone_a.new"))
  unlink(file.path(cp_dir, "clone_a.new", "manifest.tsv"))
  expect_null(dowser:::readCheckpoint(cp_dir, "a"))
})

test_that("iterateCloneDots calls functions and cuts the trees list", {
  dots <- list(mcmc_length=100,
    UCLD_SIGMA_INIT=function(cs) paste0("sigma_", cs$clone_id),
    trees=list("tree_a", "tree_b", "tree_c"))
  cs <- dplyr::tibble(clone_id="b")
  out <- dowser:::iterateCloneDots(dots, cs, c("a", "b", "c"))
  expect_equal(out$mcmc_length, 100)
  expect_equal(out$UCLD_SIGMA_INIT, "sigma_b")
  expect_equal(out$trees, list("tree_b"))
})

test_that("iterateCloneDots leaves dots without trees unchanged", {
  out <- dowser:::iterateCloneDots(list(mcmc_length=100), dplyr::tibble(clone_id="a"), "a")
  expect_equal(out, list(mcmc_length=100))
})

test_that("after an error no new clone starts, and earlier clones are saved", {
  ess_fn <- function(cid, n) if(cid == "b") NA else 1000
  s <- iterate_setup(c("a", "b", "c"), ess_fn)
  expect_error(suppressMessages(iterate_run(s)), "Clone b failed: fake failure")
  cp <- file.path(s$dir, "checkpoint_t1")
  expect_true(file.exists(file.path(cp, "clone_a", "manifest.tsv")))
  expect_false(file.exists(file.path(cp, "clone_c", "manifest.tsv")))
  expect_false(file.exists(file.path(s$dir, "t1_c.log")))
  status <- readLines(file.path(cp, "status.tsv"))
  expect_equal(status, c("clone_id\tstatus\trounds_done", "a\tdone\t1",
    "b\tfailed\t0", "c\tqueued\t0"))
})

test_that("the error message says how to resume", {
  s <- iterate_setup(c("a", "b"), function(cid, n) if(cid == "b") NA else 1000)
  expect_error(suppressMessages(iterate_run(s)), "continue=TRUE")
})

test_that("an analysis that stopped after an error can be continued", {
  fail <- TRUE
  ess_fn <- function(cid, n) if(cid == "b" && fail) NA else 1000
  s <- iterate_setup(c("a", "b", "c"), ess_fn)
  expect_error(suppressMessages(iterate_run(s)), "Clone b failed")
  fail <- FALSE
  s$calls$log <- list()
  out <- iterate_run(s, continue=TRUE)
  expect_equal(call_clones(s), c("b", "c"))
  expect_equal(out$clone_id, c("a", "b", "c"))
})

test_that("after an error in a parallel run, running clones finish and are saved, and the message says that they are awaited", {
  skip_on_os("windows")
  skip_on_cran()
  ess_fn <- function(cid, n){
    if(cid == "a"){
      Sys.sleep(3)
      1000
    }else if(cid == "b"){
      NA
    }else{
      1000
    }
  }
  s <- iterate_setup(c("a", "b", "c"), ess_fn)
  msgs <- character(0)
  expect_error(withCallingHandlers(iterate_run(s, nproc=2), message=function(m){
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }), "Clone b failed")
  cp <- file.path(s$dir, "checkpoint_t1")
  # a was running when b failed, so a was allowed to finish and was saved
  expect_true(file.exists(file.path(cp, "clone_a", "manifest.tsv")))
  # c was waiting in the queue and never started
  expect_false(file.exists(file.path(s$dir, "t1_c.log")))
  expect_true(any(grepl("No new clone will start. Waiting for 1 running clone", msgs)))
})

test_that("iterateSlimError keeps the message and drops the call", {
  e <- tryCatch(stop("some problem"), error=dowser:::iterateSlimError)
  expect_s3_class(e, "error")
  expect_equal(conditionMessage(e), "some problem")
  expect_null(conditionCall(e))
})

test_that("formatSeconds writes minutes and seconds", {
  expect_equal(dowser:::formatSeconds(c(0, 0.2, 59.7, 61, 250, NA)),
    c("-", "-", "1m00s", "1m01s", "4m10s", "-"))
})

test_that("iterateSummaryTable has one row for each clone", {
  state <- list(a=dplyr::tibble(clone_id="a", parameters=list(fake_params(100))),
    b=dplyr::tibble(clone_id="b"))
  tbl <- expect_no_warning(dowser:::iterateSummaryTable(c("a", "b"),
    c(a="unconverged", b="queued"), c(a=2L, b=0L), state, "traitfrequencies",
    200, c(a=125, b=0)))
  expect_equal(names(tbl), c("clone_id", "status", "rounds", "below_ESS", "time"))
  expect_equal(tbl$clone_id, c("a", "b"))
  expect_equal(tbl$status, c("unconverged", "queued"))
  expect_equal(tbl$rounds, c(2L, 0L))
  expect_equal(tbl$below_ESS, c(2L, NA))
  expect_equal(tbl$time, c("2m05s", "-"))
})

test_that("iterateSummaryTable gives NA for a clone with a bad parameter table", {
  state <- list(a=dplyr::tibble(clone_id="a", parameters=list(data.frame(x=1))))
  tbl <- dowser:::iterateSummaryTable("a", c(a="done"), c(a=1L), state,
    "traitfrequencies", 200, c(a=1))
  expect_true(is.na(tbl$below_ESS))
})

test_that("engine prints the summary table at the end, unless quiet", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  txt <- capture.output(iterate_run(s, quiet=0))
  expect_true(any(grepl("=== summary", txt)))
  expect_true(any(grepl("clone_id.*status.*rounds.*below_ESS.*time", txt)))
  expect_true(any(grepl("a\\s+done\\s+2\\s+0", txt)))
  expect_true(any(grepl("b\\s+done\\s+2\\s+0", txt)))
  expect_true(any(grepl("time in this call:", txt)))

  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  expect_silent(iterate_run(s, quiet=1))
})

test_that("engine prints the summary table also if the run stops with an error", {
  ess_fn <- function(cid, n) if(cid == "b") NA else 1000
  s <- iterate_setup(c("a", "b", "c"), ess_fn)
  txt <- capture.output(try(suppressMessages(iterate_run(s, quiet=0)), silent=TRUE))
  expect_true(any(grepl("a\\s+done\\s+1\\s+0", txt)))
  expect_true(any(grepl("b\\s+failed\\s+0", txt)))
  expect_true(any(grepl("c\\s+queued\\s+0", txt)))
})

test_that("the summary table shows the unconverged status and below_ESS", {
  s <- iterate_setup("a", function(cid, n) 1)
  txt <- capture.output(suppressWarnings(iterate_run(s, quiet=0, iterations=2)))
  expect_true(any(grepl("a\\s+unconverged\\s+2\\s+2", txt)))
})

test_that("parallel runs add the describe text to the start line, sequential runs do not", {
  skip_on_os("windows")
  skip_on_cran()
  describe <- function(cs) paste("info", cs$clone_id)
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  txt <- capture.output(iterate_run(s, quiet=0, nproc=2, describe=describe))
  expect_true(any(grepl("Starting clone a round 1 \\(info a\\)", txt)))
  expect_true(any(grepl("Starting clone b round 1 \\(info b\\)", txt)))
  # the text is shown only the first time that a clone starts
  expect_true(any(grepl("Starting clone a round 2", txt)))
  expect_false(any(grepl("round 2 \\(info", txt)))
  expect_equal(sum(grepl("\\(info a\\)", txt)), 1)
  expect_equal(sum(grepl("\\(info b\\)", txt)), 1)

  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  txt <- capture.output(iterate_run(s, quiet=0, describe=describe))
  expect_true(any(grepl("Starting clone a round 1", txt)))
  expect_false(any(grepl("info", txt)))
})

test_that("the describe text is shown again for a clone that starts in a later call", {
  skip_on_os("windows")
  skip_on_cran()
  describe <- function(cs) paste("info", cs$clone_id)
  s <- iterate_setup("a", function(cid, n) 1)
  suppressWarnings(capture.output(iterate_run(s, quiet=0, nproc=2, iterations=1,
    describe=describe)))
  txt <- capture.output(suppressWarnings(iterate_run(s, quiet=0, nproc=2,
    iterations=2, continue=TRUE, describe=describe)))
  expect_true(any(grepl("Starting clone a round 2 \\(info a\\)", txt)))
})

test_that("a describe function that fails does not stop the run", {
  skip_on_os("windows")
  skip_on_cran()
  s <- iterate_setup(c("a", "b"), function(cid, n) 1000)
  txt <- capture.output(out <- iterate_run(s, quiet=0, nproc=2,
    describe=function(cs) stop("no text")))
  expect_equal(out$below_ESS, c(0L, 0L))
  expect_true(any(grepl("Starting clone a round 1", txt)))
})

test_that("getEmpNucleotideFreqs gives the frequencies of A, C, G and T", {
  freqs <- dowser:::getEmpNucleotideFreqs(c("ACGT", "AACC"))
  expect_equal(as.numeric(freqs), c(0.375, 0.375, 0.125, 0.125))
  expect_equal(names(freqs), c("A", "C", "G", "T"))
})

test_that("continue uses the data of this call for a clone that had finished", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 1)
  suppressWarnings(iterate_run(s, iterations=1))
  s$clones$data <- list("a_new", "b_new")
  s$calls$log <- list()
  suppressWarnings(iterate_run(s, iterations=2, continue=TRUE))
  expect_equal(vapply(calls_of(s), function(x)x$data, character(1)), c("a_new", "b_new"))
})

test_that("continue uses the data of this call for a clone with a full checkpoint, and for a clone that did not run", {
  fail <- TRUE
  # clone b fails in the first call, and converges in its first round in the second call
  ess_fn <- function(cid, n){
    if(cid == "b") (if(fail) NA else 1000) else 1
  }
  s <- iterate_setup(c("a", "b"), ess_fn)
  suppressWarnings(iterate_run(s, iterations=1, skip_errors=TRUE))
  fail <- FALSE
  s$clones$data <- list("a_new", "b_new")
  s$calls$log <- list()
  suppressWarnings(iterate_run(s, iterations=2, continue=TRUE))
  expect_equal(vapply(calls_of(s), function(x)x$clone, character(1)), c("a", "b"))
  expect_equal(vapply(calls_of(s), function(x)x$data, character(1)), c("a_new", "b_new"))
})

test_that("the result of continue has the data of this call", {
  s <- iterate_setup("a", function(cid, n) 1)
  suppressWarnings(iterate_run(s, iterations=1))
  s$clones$data <- list("a_new")
  out <- suppressWarnings(iterate_run(s, iterations=2, continue=TRUE))
  expect_equal(out$data[[1]], "a_new")
})

# ---- the options that are not in the XML ------------------------------------------

manifest_lines <- function(s, clone) readLines(file.path(s$dir, "checkpoint_t1",
  paste0("clone_", clone), "manifest.tsv"))

test_that("the options text of a clone is saved in its checkpoint", {
  opts <- function(cs) paste0("seed=", cs$clone_id)
  ess_fn <- function(cid, n) if(cid == "b") NA else 100 * n
  s <- iterate_setup(c("a", "b"), ess_fn)
  # clone b fails, but the checkpoint of clone a is saved
  expect_error(suppressMessages(iterate_run(s, beast_options=opts)), "Clone b failed")
  expect_equal(manifest_lines(s, "a")[2], "options\tseed=a")

  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s, beast_options=opts)
  expect_equal(manifest_lines(s, "a")[2], "options\tseed=a")
  expect_equal(manifest_lines(s, "b")[2], "options\tseed=b")
})

test_that("the options text is empty if the engine gets no beast_options", {
  s <- iterate_setup("a", function(cid, n) 1000)
  iterate_run(s)
  expect_equal(manifest_lines(s, "a")[2], "options\t")
})

test_that("a clone that runs again saves the options of this call", {
  s <- iterate_setup("a", function(cid, n) 1)
  suppressWarnings(iterate_run(s, iterations=1, beast_options=function(cs) "seed=1"))
  expect_equal(manifest_lines(s, "a")[2], "options\tseed=1")
  suppressWarnings(iterate_run(s, iterations=2, continue=TRUE,
    beast_options=function(cs) "seed=2"))
  expect_equal(manifest_lines(s, "a")[2], "options\tseed=2")
})

# ==== what continue=TRUE can use ============================================================

# ---- getClonesToRestart ------------------------------------------------------------

xml_md5_of <- function(text){
  path <- tempfile(fileext=".xml")
  writeLines(text, path)
  unname(tools::md5sum(path))
}

# An analysis with a checkpoint for each clone. xml_text gives the XML text of each clone.
cp_setup <- function(ids, xml_text, options="seed=1"){
  s <- iterate_setup(ids, function(cid, n) 1000)
  for(cid in ids){
    if(!is.null(xml_text[[cid]])){
      writeLines(xml_text[[cid]], file.path(s$dir, paste0("t1_", cid, ".xml")))
    }
    result <- s$work(s$clones[match(cid, ids), ], NULL)
    dowser:::writeIterateCheckpoint(s$dir, "t1", ids, cid, 1L, result, options)
  }
  s
}
same_options <- function(cid) "seed=1"
never <- function(cid) stop("must not be called")

test_that("getClonesToRestart has nothing to do if there is no checkpoint folder", {
  out <- dowser:::getClonesToRestart(tempfile("no_dir_"), "t1", c("a", "b"), never, never)
  expect_equal(out, character(0))
})

test_that("getClonesToRestart accepts the same options and the same XML, without a message", {
  s <- cp_setup(c("a", "b"), list(a="<x a/>", b="<x b/>"))
  md5 <- c(a=xml_md5_of("<x a/>"), b=xml_md5_of("<x b/>"))
  out <- expect_no_message(dowser:::getClonesToRestart(s$dir, "t1", c("a", "b"),
    same_options, function(cid) md5[[cid]]))
  expect_equal(out, character(0))
})

test_that("getClonesToRestart restarts a clone with other options, without making its XML", {
  s <- cp_setup(c("a", "b"), list(a="<x a/>", b="<x b/>"))
  expect_message(
    out <- dowser:::getClonesToRestart(s$dir, "t1", c("a", "b"),
      function(cid) "seed=2", never),
    "clone(s) a, b changed since the earlier run", fixed=TRUE)
  expect_equal(out, c("a", "b"))
})

test_that("getClonesToRestart compares the options of each clone", {
  s <- cp_setup(c("a", "b"), list(a="<x a/>", b="<x b/>"))
  md5 <- c(a=xml_md5_of("<x a/>"), b=xml_md5_of("<x b/>"))
  out <- suppressMessages(dowser:::getClonesToRestart(s$dir, "t1", c("a", "b"),
    function(cid) if(cid == "b") "seed=2" else "seed=1", function(cid) md5[[cid]]))
  expect_equal(out, "b")
})

test_that("getClonesToRestart restarts only the clone with another XML", {
  s <- cp_setup(c("a", "b", "c"), list(a="<x a/>", b="<x b/>", c="<x c/>"))
  md5 <- c(a=xml_md5_of("<x a/>"), b=xml_md5_of("<x b changed/>"), c=xml_md5_of("<x c/>"))
  expect_message(
    out <- dowser:::getClonesToRestart(s$dir, "t1", c("a", "b", "c"),
      same_options, function(cid) md5[[cid]]),
    "clone(s) b changed since the earlier run", fixed=TRUE)
  expect_equal(out, "b")
})

test_that("the message says what is deleted and that the other clones continue", {
  s <- cp_setup(c("a", "b"), list(a="<x a/>", b="<x b/>"))
  msgs <- character(0)
  withCallingHandlers(
    dowser:::getClonesToRestart(s$dir, "t1", c("a", "b"), same_options,
      function(cid) xml_md5_of("<other/>")),
    message=function(m){ msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("checkpoints and BEAST files are deleted", msgs, fixed=TRUE)))
  expect_true(any(grepl("start again from round 1", msgs, fixed=TRUE)))
  expect_true(any(grepl("The other clones continue", msgs, fixed=TRUE)))
})

test_that("getClonesToRestart does not restart a clone if the XML cannot be made", {
  s <- cp_setup(c("a", "b"), list(a="<x a/>", b="<x b/>"))
  out <- dowser:::getClonesToRestart(s$dir, "t1", c("a", "b"), same_options,
    function(cid) NA_character_)
  expect_equal(out, character(0))
})

test_that("getClonesToRestart does not compare a clone without checkpoint or without XML", {
  s <- cp_setup(c("a", "b"), list(a="<x a/>"))   # no XML file for clone b
  asked <- character(0)
  out <- dowser:::getClonesToRestart(s$dir, "t1", c("a", "b", "c"), same_options,
    function(cid){
      asked <<- c(asked, cid)
      xml_md5_of("<x a/>")
    })
  expect_equal(asked, "a")   # b has no XML in the manifest, c has no checkpoint
  expect_equal(out, character(0))
})

# the format of an earlier version: a state line, and no kind and size for each file
old_manifest <- function(s){
  manifest <- file.path(s$dir, "checkpoint_t1", "clone_a", "manifest.tsv")
  writeLines(c("round\t1", "state\tfull", "options\tseed=1", "t1_a.xml\tabc", "result.rds\tdef"),
    manifest)
}

test_that("getClonesToRestart does not compare a checkpoint with an old manifest", {
  s <- cp_setup("a", list(a="<x a/>"))
  old_manifest(s)
  out <- dowser:::getClonesToRestart(s$dir, "t1", "a", never, never)
  expect_equal(out, character(0))
})

test_that("a clone with an old manifest is reported as damaged by continue, and starts again", {
  s <- cp_setup("a", list(a="<x a/>"))
  old_manifest(s)
  s$calls$log <- list()
  expect_message(suppressWarnings(iterate_run(s, iterations=1, continue=TRUE)), "is damaged")
  expect_null(calls_of(s)[[1]]$resume)
})

# ============================================================================================
# Part 3: getTimeTreesIterate, backup and restoreIterateBackup, readFromLastCheckpoint
# ============================================================================================

make_airr_clones <- function(ids, trait_values, sequences=NULL){
  data <- lapply(seq_along(ids), function(i){
    methods::new("airrClone",
      data=data.frame(sequence_id=paste0("s", seq_along(trait_values[[i]])),
        celltype=trait_values[[i]],
        sequence=if(is.null(sequences)){
          rep(c("ACGT", "AACC"), length.out=length(trait_values[[i]]))
        }else{
          sequences[[i]]
        },
        stringsAsFactors=FALSE),
      clone=ids[i])
  })
  dplyr::tibble(clone_id=ids, data=data)
}

# Fake getTimeTrees. It writes an XML file that depends on the data of the clone
# and on the arguments that define the chain. If it gets the argument xml_check,
# it only writes the XML and calls xml_check (this is how getTimeTreesIterate
# gets the XML without running BEAST). Otherwise it adds a line to the .log
# file of the clone. The number of lines is the round. All arguments of the
# call are recorded.
fake_getTimeTrees <- function(calls, ess_fn){
  function(clones, resume_clones=NULL, dir, id, ...){
    args <- list(...)
    cid <- as.character(clones$clone_id)
    # the seed, the java flag and the BEAST path are on the BEAST command line, not in the XML
    xml_args <- args[setdiff(names(args), c("xml_check", "quiet", "nproc", "burnin",
      "asr", "trim_ids", "low_ram", "posterior", "template", "beast", "time", "seed", "java"))]
    template_text <- if(is.character(args$template) && file.exists(args$template)){
      readLines(args$template)
    }else{
      ""
    }
    xml <- file.path(dir, paste0(id, "_", cid, ".xml"))
    writeLines(c("<beast>", template_text, paste(clones$data[[1]]@data$sequence, collapse=","),
      paste(deparse(xml_args), collapse=" "), "</beast>"), xml)
    if(!is.null(args$xml_check)){
      args$xml_check(xml)
      return(clones)
    }
    f <- file.path(dir, paste0(id, "_", cid, ".log"))
    cat("round\n", file=f, append=TRUE)
    n <- length(readLines(f))
    calls$log[[length(calls$log) + 1]] <- c(
      list(clone=cid, round=n, resume=resume_clones, nrow=nrow(clones),
        dir=dir, id=id), args)
    ess <- ess_fn(cid, n)
    if(is.na(ess)){
      stop("fake failure")
    }
    clones$parameters <- list(fake_params(ess))
    clones
  }
}

# run f and return the value and all messages
with_messages <- function(f){
  msgs <- character(0)
  value <- withCallingHandlers(f(), message=function(m){
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  list(value=value, messages=msgs)
}

run_wrapper <- function(clones, dir, ...){
  dowser::getTimeTreesIterate(clones, template="template.xml", beast="beast_bin",
    dir=dir, id="t1", time="time", quiet=1, ...)
}

new_calls <- function(){
  calls <- new.env()
  calls$log <- list()
  calls
}

test_that("every clone gets the same trait_list, without ?", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list(c("GC", "?"), c("GC", "MBC")))
  run_wrapper(clones, tempfile("wrapper_test_"), trait="celltype")
  expect_equal(length(calls$log), 2)
  for(x in calls$log){
    expect_equal(x$trait_list, c("GC", "MBC"))
    expect_equal(x$trait, "celltype")
  }
})

test_that("a trait_list from the user is passed on without change", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones("a", list(c("GC", "MBC")))
  run_wrapper(clones, tempfile("wrapper_test_"), trait="celltype",
    trait_list=c("MBC", "GC", "PC"))
  expect_equal(calls$log[[1]]$trait_list, c("MBC", "GC", "PC"))
})

test_that("there is no trait_list if there is no trait", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones("a", list(c("GC", "MBC")))
  run_wrapper(clones, tempfile("wrapper_test_"))
  expect_null(calls$log[[1]]$trait_list)
})

test_that("each call has one clone, posterior none, one core and the first call is not a resume", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 100 * n), .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  out <- run_wrapper(clones, tempfile("wrapper_test_"))
  expect_equal(length(calls$log), 4)
  for(x in calls$log){
    expect_equal(x$nrow, 1)
    expect_equal(x$posterior, "none")
    expect_equal(x$nproc, 1)
    expect_equal(x$template, "template.xml")
    expect_equal(x$time, "time")
    if(x$round == 1) expect_null(x$resume) else expect_equal(x$resume, x$clone)
  }
  expect_equal(out$below_ESS, c(0L, 0L))
  expect_equal(out$clone_id, c("a", "b"))
})

test_that("arguments in ... are passed on, and functions are called for each clone", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  run_wrapper(clones, tempfile("wrapper_test_"), mcmc_length=500,
    UCLD_SIGMA_INIT=function(cs) paste0("sigma_", cs$clone_id))
  by_clone <- stats::setNames(calls$log, vapply(calls$log, function(x)x$clone, character(1)))
  expect_equal(by_clone$a$UCLD_SIGMA_INIT, "sigma_a")
  expect_equal(by_clone$b$UCLD_SIGMA_INIT, "sigma_b")
  expect_equal(by_clone$a$mcmc_length, 500)
})

test_that("starting trees are cut to the tree of each clone", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  run_wrapper(clones, tempfile("wrapper_test_"), trees=list("tree_a", "tree_b"))
  by_clone <- stats::setNames(calls$log, vapply(calls$log, function(x)x$clone, character(1)))
  expect_equal(by_clone$a$trees, list("tree_a"))
  expect_equal(by_clone$b$trees, list("tree_b"))
})

test_that("readBEAST is not called when posterior is none", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  rb <- new.env()
  rb$n <- 0
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000),
    readBEAST=function(clones, ...){ rb$n <- rb$n + 1; clones },
    .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  run_wrapper(clones, tempfile("wrapper_test_"))
  expect_equal(rb$n, 0)
})

test_that("readBEAST is called once at the end when a posterior is wanted, and below_ESS is made again", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  rb <- new.env()
  rb$n <- 0
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000),
    readBEAST=function(clones, posterior, ...){
      rb$n <- rb$n + 1
      rb$posterior <- posterior
      rb$nrow <- nrow(clones)
      clones$parameters <- lapply(seq_len(nrow(clones)), function(i) fake_params(1))
      clones
    },
    .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  out <- run_wrapper(clones, tempfile("wrapper_test_"), posterior="all")
  expect_equal(rb$n, 1)
  expect_equal(rb$posterior, "all")
  expect_equal(rb$nrow, 2)
  # the fake readBEAST gives ESS 1, so both parameters are below the cutoff
  expect_equal(out$below_ESS, c(2L, 2L))
})

test_that("dir must be given", {
  clones <- make_airr_clones("a", list("GC"))
  expect_error(dowser::getTimeTreesIterate(clones, template="t.xml", beast="b",
    dir=NULL, id="t1", time="time"), "dir must be specified")
})

test_that("there is no call when there are no clones, and below_ESS is in the result", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones(character(0), list())
  out <- run_wrapper(clones, tempfile("wrapper_test_"))
  expect_equal(length(calls$log), 0)
  expect_true("below_ESS" %in% names(out))
})

test_that("backup=TRUE moves the old files, and the backup stays after the run", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  run_wrapper(clones, dir)
  expect_equal(length(readLines(file.path(dir, "t1_a.log"))), 1)
  # a second run without continue: the old file has 1 line. A moved file is
  # not there, so the new file also has 1 line, and not 2.
  suppressMessages(run_wrapper(clones, dir, backup=TRUE))
  expect_equal(length(readLines(file.path(dir, "t1_a.log"))), 1)
  expect_equal(length(list.files(dir, pattern="^backup_t1_")), 1)
})

test_that("the backup can be put back, and it is deleted when it is used", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 100 * n), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  run_wrapper(clones, dir)
  old_log <- readLines(file.path(dir, "t1_a.log"))
  expect_equal(length(old_log), 2)
  suppressMessages(run_wrapper(clones, dir, backup=TRUE))
  expect_equal(length(list.files(dir, pattern="^backup_t1_")), 1)
  # the log of the new run is not the same file as the old one
  writeLines("new", file.path(dir, "t1_a.log"))
  suppressMessages(restoreIterateBackup(dir, "t1"))
  expect_equal(readLines(file.path(dir, "t1_a.log")), old_log)
  expect_equal(length(list.files(dir, pattern="^backup_t1_")), 0)
})

test_that("backup=TRUE with continue=TRUE copies the files, the run goes on, and the backup stays", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  run_wrapper(clones, dir)
  old_log <- readLines(file.path(dir, "t1_a.log"))
  res <- with_messages(function() run_wrapper(clones, dir, backup=TRUE, continue=TRUE))
  expect_true(any(grepl("copied to", res$messages)))
  # the files are still in dir, and the backup has a copy
  expect_equal(readLines(file.path(dir, "t1_a.log")), old_log)
  bk <- list.files(dir, pattern="^backup_t1_", full.names=TRUE)
  expect_equal(length(bk), 1)
  expect_equal(readLines(file.path(bk, "t1_a.log")), old_log)
  # every call with backup=TRUE makes a new backup folder
  suppressMessages(run_wrapper(clones, dir, backup=TRUE, continue=TRUE))
  expect_equal(length(list.files(dir, pattern="^backup_t1_")), 2)
})

test_that("the backup is kept if the run stops with an error", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  run_wrapper(clones, dir)
  testthat::local_mocked_bindings(
    getTimeTrees=function(...) stop("BEAST failed"), .package="dowser")
  expect_error(suppressMessages(run_wrapper(clones, dir, backup=TRUE)), "BEAST failed")
  expect_equal(length(list.files(dir, pattern="^backup_t1_")), 1)
})

test_that("a parallel run shows the nucleotide frequencies if the template uses them", {
  skip_if_not_installed("testthat", "3.2.0")
  skip_on_os("windows")
  skip_on_cran()
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1000), .package="dowser")
  clones <- make_airr_clones("a", list(c("GC", "MBC")))
  with_freqs <- tempfile(fileext=".xml")
  writeLines(c("<beast>", "<frequencies>${EMP_EQFREQS}</frequencies>", "</beast>"), with_freqs)
  without_freqs <- tempfile(fileext=".xml")
  writeLines(c("<beast>", "</beast>"), without_freqs)
  run <- function(template){
    capture.output(dowser::getTimeTreesIterate(clones, template=template,
      beast="beast_bin", dir=tempfile("wrapper_test_"), id="t1", time="time",
      quiet=0, nproc=2))
  }
  txt <- run(with_freqs)
  expect_true(any(grepl("nucleotide frequencies A C G T: 0.375 0.375 0.125 0.125", txt)))
  txt <- run(without_freqs)
  expect_true(any(grepl("Starting clone a round 1", txt)))
  expect_false(any(grepl("nucleotide frequencies", txt)))
})

test_that("the nucleotide frequencies are not rounded, and are shown only once for each clone", {
  skip_if_not_installed("testthat", "3.2.0")
  skip_on_os("windows")
  skip_on_cran()
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 100 * n), .package="dowser")
  # A 2, C 1, G 1, T 2 gives the frequencies 1/3, 1/6, 1/6, 1/3
  clones <- make_airr_clones("a", list(c("GC", "MBC")), sequences=list(c("AAC", "GTT")))
  tmpl <- tempfile(fileext=".xml")
  writeLines(c("<beast>", "<f>${EMP_EQFREQS}</f>", "</beast>"), tmpl)
  txt <- capture.output(dowser::getTimeTreesIterate(clones, template=tmpl,
    beast="beast_bin", dir=tempfile("wrapper_test_"), id="t1", time="time",
    quiet=0, nproc=2))
  expect_equal(sum(grepl("nucleotide frequencies", txt)), 1)
  expect_true(any(grepl(paste("A C G T:", 2/6, 1/6, 1/6, 2/6), txt, fixed=TRUE)))
  expect_true(any(grepl("Starting clone a round 2", txt)))
})

# ---- continue only uses checkpoints that match the call -----------------------------

test_that("arguments that do not change the chain can change on continue", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=1))
  calls$log <- list()
  suppressWarnings(run_wrapper(clones, dir, iterations=3, continue=TRUE, burnin=20,
    ess_cutoff=500, ignore=c("freqParameter"), skip_errors=TRUE, nproc=1))
  expect_equal(length(calls$log), 2)
})

test_that("continue with the same arguments continues all clones", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, mcmc_length=1000, seed=7, iterations=1))
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir,
    mcmc_length=1000, seed=7, iterations=2, continue=TRUE)))
  expect_equal(length(res$messages), 0)
  expect_equal(length(calls$log), 2)
  for(x in calls$log){
    expect_equal(x$resume, x$clone)
    expect_equal(x$round, 2)
  }
})

test_that("pre-flight XML calls are not counted as rounds and write no file in dir", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=1))
  xml_before <- tools::md5sum(file.path(dir, "t1_a.xml"))
  calls$log <- list()
  suppressWarnings(run_wrapper(clones, dir, iterations=1, continue=TRUE))
  expect_equal(length(calls$log), 0)
  expect_equal(tools::md5sum(file.path(dir, "t1_a.xml")), xml_before)
  expect_equal(length(list.files(tempdir(), pattern="^iterate_xml_")), 0)
})

# A runner with another BEAST path
run_beast <- function(clones, dir, beast, ...){
  dowser::getTimeTreesIterate(clones, template="template.xml", beast=beast, dir=dir,
    id="t1", time="time", quiet=1, ...)
}
restart_message <- function(ids) paste0("clone(s) ", paste(ids, collapse=", "),
  " changed since the earlier run")

test_that("a changed seed, java flag or BEAST path starts the clones again from round 1", {
  skip_if_not_installed("testthat", "3.2.0")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  changes <- list(
    list(change=list(seed=8), saved="seed=8 "),
    list(change=list(java=FALSE), saved="java=FALSE"),
    list(change=list(beast="other_beast"), saved="other_beast$"))
  for(ch in changes){
    calls <- new_calls()
    testthat::local_mocked_bindings(
      getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
    dir <- tempfile("wrapper_test_")
    suppressWarnings(run_beast(clones, dir, "beast_bin", iterations=1, seed=7))
    calls$log <- list()
    args <- list(clones=clones, dir=dir, beast="beast_bin", iterations=1, continue=TRUE, seed=7)
    for(n in names(ch$change)) args[[n]] <- ch$change[[n]]
    res <- with_messages(function() suppressWarnings(do.call(run_beast, args)))
    info <- names(ch$change)
    expect_true(any(grepl(restart_message(c("a", "b")), res$messages, fixed=TRUE)), info=info)
    expect_equal(length(calls$log), 2, info=info)
    for(x in calls$log) expect_null(x$resume, info=info)
    manifest <- readLines(file.path(dir, "checkpoint_t1", "clone_a", "manifest.tsv"))
    expect_match(manifest[2], ch$saved, info=info)
  }
})

test_that("a seed that is removed is a changed option", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=1, seed=7))
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir, iterations=1,
    continue=TRUE)))
  expect_true(any(grepl(restart_message("a"), res$messages, fixed=TRUE)))
  expect_null(calls$log[[1]]$resume)
})

test_that("a changed argument that is in the XML starts all clones again from round 1", {
  skip_if_not_installed("testthat", "3.2.0")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  for(change in list(list(mcmc_length=2000), list(KAPPA_PRIOR_M=0.9), list(log_target=200))){
    calls <- new_calls()
    testthat::local_mocked_bindings(
      getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
    dir <- tempfile("wrapper_test_")
    base <- list(mcmc_length=1000, KAPPA_PRIOR_M=0.67, log_target=100)
    suppressWarnings(do.call(run_wrapper, c(list(clones, dir, iterations=1), base)))
    calls$log <- list()
    args <- base
    for(n in names(change)) args[[n]] <- change[[n]]
    res <- with_messages(function() suppressWarnings(do.call(run_wrapper,
      c(list(clones, dir, iterations=1, continue=TRUE), args))))
    expect_true(any(grepl(restart_message(c("a", "b")), res$messages, fixed=TRUE)),
      info=names(change))
    expect_equal(length(calls$log), 2, info=names(change))
    for(x in calls$log) expect_null(x$resume, info=names(change))
  }
})

test_that("a changed template file starts all clones again from round 1", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  tmpl <- tempfile(fileext=".xml")
  writeLines("<beast a='1'/>", tmpl)
  run <- function(...){
    dowser::getTimeTreesIterate(clones, template=tmpl, beast="beast_bin", dir=dir,
      id="t1", time="time", quiet=1, ...)
  }
  suppressWarnings(run(iterations=1))
  writeLines("<beast a='2'/>", tmpl)
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run(iterations=1, continue=TRUE)))
  expect_true(any(grepl(restart_message("a"), res$messages, fixed=TRUE)))
  expect_null(calls$log[[1]]$resume)
})

test_that("a function argument with the same values continues, and with other values restarts the clones", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=1, SIGMA=function(cs) 1))
  # other code, same values: the XML is the same
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir, iterations=2,
    continue=TRUE, SIGMA=function(cs) {1})))
  expect_equal(length(res$messages), 0)
  for(x in calls$log) expect_equal(x$resume, x$clone)
  # other values: the XML is not the same, so the clones start again
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir, iterations=1,
    continue=TRUE, SIGMA=function(cs) 2)))
  expect_true(any(grepl(restart_message(c("a", "b")), res$messages, fixed=TRUE)))
  for(x in calls$log) expect_null(x$resume)
})

test_that("a function and a number with the same value give the same XML, so the clones continue", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=1, SIGMA=function(cs) 0.5))
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir, iterations=2,
    continue=TRUE, SIGMA=0.5)))
  expect_equal(length(res$messages), 0)
  expect_equal(calls$log[[1]]$resume, "a")
})

test_that("changed data of one clone restarts only that clone", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  dir <- tempfile("wrapper_test_")
  clones1 <- make_airr_clones(c("a", "b"), list("GC", "GC"), sequences=list("ACGT", "ACGT"))
  suppressWarnings(run_wrapper(clones1, dir, iterations=1))
  clones2 <- make_airr_clones(c("a", "b"), list("GC", "GC"), sequences=list("ACGT", "TTTT"))
  # files of another analysis and a file with another name must not be deleted
  file.create(file.path(dir, c("t2_b.log", "notes.txt")))
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones2, dir, iterations=2,
    continue=TRUE)))
  expect_true(all(file.exists(file.path(dir, c("t2_b.log", "notes.txt")))))
  expect_true(any(grepl(restart_message("b"), res$messages, fixed=TRUE)))
  expect_false(any(grepl("clone(s) a", res$messages, fixed=TRUE)))
  a_call <- Filter(function(x)x$clone == "a", calls$log)[[1]]
  b_call <- Filter(function(x)x$clone == "b", calls$log)[[1]]
  # a continues with round 2
  expect_equal(a_call$resume, "a")
  expect_equal(a_call$round, 2)
  # b starts again: no resume, and its old files were deleted, so this is the first line of its log
  expect_null(b_call$resume)
  expect_equal(b_call$round, 1)
  # b has a new checkpoint, and the checkpoint folder has no leftover
  expect_true(file.exists(file.path(dir, "checkpoint_t1", "clone_b", "manifest.tsv")))
  expect_equal(length(list.files(file.path(dir, "checkpoint_t1"), pattern="\\.new$")), 0)
})

test_that("a failed clone and a changed seed: the clone with a checkpoint starts again, and the failed clone starts", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  fail <- TRUE
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) if(cid == "b" && fail) NA else 1),
    .package="dowser")
  clones <- make_airr_clones(c("a", "b"), list("GC", "GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=1, skip_errors=TRUE, seed=7))
  fail <- FALSE
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir, iterations=1,
    skip_errors=TRUE, continue=TRUE, seed=8)))
  expect_true(any(grepl(restart_message("a"), res$messages, fixed=TRUE)))
  expect_equal(vapply(calls$log, function(x)x$clone, ""), c("a", "b"))
  for(x in calls$log) expect_null(x$resume)
})

test_that("with backup=TRUE, a continued run that restarts a clone keeps a copy of the old files", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=2, seed=7))
  old_log <- readLines(file.path(dir, "t1_a.log"))
  expect_equal(length(old_log), 2)
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir, iterations=1,
    continue=TRUE, seed=8, backup=TRUE)))
  expect_true(any(grepl("copied to", res$messages)))
  bk <- list.files(dir, pattern="^backup_t1_", full.names=TRUE)
  expect_equal(length(bk), 1)    # the backup has the old files of the clone that started again
  expect_equal(readLines(file.path(bk, "t1_a.log")), old_log)
  manifest <- readLines(file.path(bk, "checkpoint_t1", "clone_a", "manifest.tsv"))
  expect_match(manifest[2], "seed=7", fixed=TRUE)
  # the clone started again with fresh files
  expect_equal(calls$log[[1]]$round, 1)
  expect_null(calls$log[[1]]$resume)
})

test_that("a continued run with backup=TRUE can be undone with restoreIterateBackup", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=2, seed=7))
  old_log <- readLines(file.path(dir, "t1_a.log"))
  suppressMessages(suppressWarnings(run_wrapper(clones, dir, iterations=1, continue=TRUE,
    seed=8, backup=TRUE)))
  suppressMessages(restoreIterateBackup(dir, "t1"))
  expect_equal(readLines(file.path(dir, "t1_a.log")), old_log)
  # the old analysis continues with its own options
  calls$log <- list()
  res <- with_messages(function() suppressWarnings(run_wrapper(clones, dir, iterations=3,
    continue=TRUE, seed=7)))
  expect_equal(length(res$messages), 0)
  expect_equal(calls$log[[1]]$resume, "a")
  expect_equal(calls$log[[1]]$round, 3)
})

test_that("a restart for another option gives no warning", {
  skip_if_not_installed("testthat", "3.2.0")
  calls <- new_calls()
  testthat::local_mocked_bindings(
    getTimeTrees=fake_getTimeTrees(calls, function(cid, n) 1), .package="dowser")
  clones <- make_airr_clones("a", list("GC"))
  dir <- tempfile("wrapper_test_")
  suppressWarnings(run_wrapper(clones, dir, iterations=1, seed=7))
  expect_no_warning(suppressMessages(run_wrapper(clones, dir, iterations=1, continue=TRUE,
    seed=8, ess_cutoff=0)))
})

# ---- the list of types of a trait (makeTraitList). "?" is an ambiguous tip value and is never
# ---- counted as a type.

test_that("makeTraitList sorts, removes duplicates and removes ?", {
  values <- list(c("GC", "?", "GC"), c("MBC", "GC", "?"))
  expect_equal(dowser:::makeTraitList(values), c("GC", "MBC"))
})

test_that("makeTraitList gives the same result when there is no ?", {
  values <- list(c("GC", "MBC"), "PC")
  expect_equal(dowser:::makeTraitList(values), c("GC", "MBC", "PC"))
})

# ==== backup and restoreIterateBackup ===========================================================

test_that("makeBackup moves the BEAST files and the checkpoint folder", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s)
  backup_dir <- dowser:::makeBackup(s$dir, "t1", c("a", "b"), FALSE)
  expect_true(dir.exists(backup_dir))
  expect_equal(list.files(s$dir), basename(backup_dir))
  moved <- list.files(backup_dir)
  expect_true(all(c("README.txt", "checkpoint_t1", "t1_a.log", "t1_a.trees",
    "t1_a.xml.state", "t1_b.log") %in% moved))
})

test_that("makeBackup does not move files of other analyses or clones", {
  s <- iterate_setup("a", function(cid, n) 1000)
  iterate_run(s)
  file.create(file.path(s$dir, c("t2_a.log", "t1_other.log", "notes.txt")))
  backup_dir <- dowser:::makeBackup(s$dir, "t1", "a", FALSE)
  expect_true(all(file.exists(file.path(s$dir, c("t2_a.log", "t1_other.log", "notes.txt")))))
  expect_false("t1_other.log" %in% list.files(backup_dir))
})

test_that("makeBackup gives NULL if there is nothing to move", {
  dir <- tempfile("backup_test_")
  dir.create(dir)
  expect_null(dowser:::makeBackup(dir, "t1", c("a", "b"), FALSE))
  expect_equal(length(list.files(dir)), 0)
})

test_that("makeBackup makes a new folder name if the name is taken", {
  s <- iterate_setup("a", function(cid, n) 1000)
  iterate_run(s)
  first <- dowser:::makeBackup(s$dir, "t1", "a", FALSE)
  file.create(file.path(s$dir, "t1_a.log"))
  second <- dowser:::makeBackup(s$dir, "t1", "a", FALSE)
  expect_false(first == second)
  expect_true(dir.exists(first))
  expect_true(dir.exists(second))
})

test_that("restoreIterateBackup puts the old files back", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s)
  old_log <- readLines(file.path(s$dir, "t1_a.log"))
  backup_dir <- dowser:::makeBackup(s$dir, "t1", c("a", "b"), FALSE)
  # files of a new run that was started by mistake
  writeLines("new run", file.path(s$dir, "t1_a.log"))
  dir.create(file.path(s$dir, "checkpoint_t1"))
  file.create(file.path(s$dir, "checkpoint_t1", "new_marker"))

  restored <- suppressMessages(restoreIterateBackup(s$dir, "t1"))
  expect_true(all(c("t1_a.log", "checkpoint_t1") %in% restored))
  expect_equal(readLines(file.path(s$dir, "t1_a.log")), old_log)
  expect_false(file.exists(file.path(s$dir, "checkpoint_t1", "new_marker")))
  expect_true(file.exists(file.path(s$dir, "checkpoint_t1", "clone_a", "manifest.tsv")))
  expect_false(dir.exists(backup_dir))
  expect_false(file.exists(file.path(s$dir, "README.txt")))
})

test_that("a restored run can be continued", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 1)
  suppressWarnings(iterate_run(s, iterations=1))
  dowser:::makeBackup(s$dir, "t1", c("a", "b"), FALSE)
  writeLines("new run", file.path(s$dir, "t1_a.log"))
  suppressMessages(restoreIterateBackup(s$dir, "t1"))
  s$calls$log <- list()
  suppressWarnings(iterate_run(s, iterations=2, continue=TRUE))
  a_call <- Filter(function(x)x$clone == "a", calls_of(s))[[1]]
  expect_equal(a_call$n, 2)
})

test_that("restoreIterateBackup uses the newest backup of the id", {
  dir <- tempfile("backup_test_")
  dir.create(dir)
  make_backup <- function(name, text){
    path <- file.path(dir, name)
    dir.create(path)
    writeLines(text, file.path(path, "t1_a.log"))
  }
  make_backup("backup_t1_20260101_000000", "old1")
  make_backup("backup_t1_20260102_000000", "old2")
  make_backup("backup_t1_x_20260103_000000", "other id")
  suppressMessages(restoreIterateBackup(dir, "t1"))
  expect_equal(readLines(file.path(dir, "t1_a.log")), "old2")
  expect_true(dir.exists(file.path(dir, "backup_t1_20260101_000000")))
  expect_true(dir.exists(file.path(dir, "backup_t1_x_20260103_000000")))
})

test_that("restoreIterateBackup can use a given backup folder", {
  dir <- tempfile("backup_test_")
  dir.create(dir)
  path <- file.path(dir, "my_backup")
  dir.create(path)
  writeLines("kept", file.path(path, "t1_a.log"))
  suppressMessages(restoreIterateBackup(dir, "t1", backup=path))
  expect_equal(readLines(file.path(dir, "t1_a.log")), "kept")
})

test_that("restoreIterateBackup stops if there is no backup", {
  dir <- tempfile("backup_test_")
  dir.create(dir)
  expect_error(restoreIterateBackup(dir, "t1"), "No backup folder found")
  expect_error(restoreIterateBackup(dir, "t1", backup=file.path(dir, "none")),
    "not found")
})


test_that("makeBackup makes the backup and says where it is", {
  s <- iterate_setup("a", function(cid, n) 1000)
  iterate_run(s)
  expect_message(res <- dowser:::makeBackup(s$dir, "t1", "a", FALSE),
    "moved to")
  expect_true(dir.exists(res))
  expect_false(file.exists(file.path(s$dir, "t1_a.log")))
})

test_that("a new run after a backup has no old checkpoint to delete", {
  s <- iterate_setup("a", function(cid, n) 1000)
  iterate_run(s)
  suppressMessages(dowser:::makeBackup(s$dir, "t1", "a", FALSE))
  res <- with_warnings(function() iterate_run(s))
  expect_false(any(grepl("Deleting the existing checkpoints", res$warnings)))
})

test_that("moveFiles can copy: the original stays", {
  dir <- tempfile("move_test_")
  dir.create(dir)
  writeLines("text", file.path(dir, "f.txt"))
  dir.create(file.path(dir, "folder"))
  writeLines("inside", file.path(dir, "folder", "g.txt"))
  dir.create(file.path(dir, "target"))
  dowser:::moveFiles(file.path(dir, "f.txt"), file.path(dir, "target", "f.txt"), keep=TRUE)
  dowser:::moveFiles(file.path(dir, "folder"), file.path(dir, "target", "folder"), keep=TRUE)
  expect_true(file.exists(file.path(dir, "f.txt")))
  expect_true(file.exists(file.path(dir, "folder", "g.txt")))
  expect_equal(readLines(file.path(dir, "target", "f.txt")), "text")
  expect_equal(readLines(file.path(dir, "target", "folder", "g.txt")), "inside")
  # without keep, the original is gone
  dowser:::moveFiles(file.path(dir, "f.txt"), file.path(dir, "moved.txt"))
  expect_false(file.exists(file.path(dir, "f.txt")))
  expect_true(file.exists(file.path(dir, "moved.txt")))
})

test_that("makeBackup with contine=TRUE copies the files and the checkpoint folder", {
  s <- iterate_setup(c("a", "b"), function(cid, n) 100 * n)
  iterate_run(s)
  before <- sort(list.files(s$dir))
  backup_dir <- dowser:::makeBackup(s$dir, "t1", c("a", "b"), TRUE)
  expect_equal(sort(list.files(s$dir)), sort(c(before, basename(backup_dir))))
  expect_true(all(c("README.txt", "checkpoint_t1", "t1_a.log", "t1_b.trees") %in%
    list.files(backup_dir)))
  expect_equal(unname(tools::md5sum(file.path(backup_dir, "t1_a.log"))),
    unname(tools::md5sum(file.path(s$dir, "t1_a.log"))))
  expect_true(file.exists(file.path(backup_dir, "checkpoint_t1", "clone_a", "manifest.tsv")))
})

# ==== readFromLastCheckpoint ================================================================

# A run that finished 2 rounds of each clone and saved the checkpoints, and then a third
# round that did not end: each file has a third line, and there is no checkpoint for it.
unfinished_run <- function(ids=c("a", "b")){
  s <- iterate_setup(ids, function(cid, n) 1000)
  for(cid in ids){
    i <- match(cid, ids)
    s$work(s$clones[i, ], NULL)
    result <- s$work(s$clones[i, ], cid)
    dowser:::writeIterateCheckpoint(s$dir, "t1", ids, cid, 2L, result, "seed=1")
    s$work(s$clones[i, ], cid)
  }
  s
}

# readBEAST that records the folder, its files and the arguments of the call
spy_readBEAST <- function(seen, fail=FALSE){
  function(clones, dir, id, beast, burnin, trait, nproc, quiet, posterior, asr, low_ram, trim_ids){
    seen$dir <- dir
    seen$files <- sort(list.files(dir))
    seen$content <- lapply(stats::setNames(seen$files, seen$files),
      function(f) readLines(file.path(dir, f), warn=FALSE))
    seen$args <- list(id=id, beast=beast, burnin=burnin, trait=trait, nproc=nproc, quiet=quiet,
      posterior=posterior, asr=asr, low_ram=low_ram, trim_ids=trim_ids,
      clone_ids=as.character(clones$clone_id))
    if(fail) stop("readBEAST failed")
    clones$trees <- as.list(rep(NA, nrow(clones)))
    clones$parameters <- lapply(seq_len(nrow(clones)), function(i) fake_params(1000))
    clones
  }
}

test_that("iterateCopyPrefix copies the first bytes of a file", {
  from <- tempfile("prefix_test_")
  writeBin(charToRaw("0123456789"), from)
  to <- tempfile("prefix_test_")
  dowser:::iterateCopyPrefix(from, to, 7, chunk=3)
  expect_equal(rawToChar(readBin(to, "raw", 100)), "0123456")
  expect_equal(file.size(from), 10)
  dowser:::iterateCopyPrefix(from, to, 0)
  expect_equal(file.size(to), 0)
  dowser:::iterateCopyPrefix(from, to, 100)    # more than the file has: the whole file
  expect_equal(file.size(to), 10)
})

test_that("the files that readBEAST needs are rebuilt as they were at the end of the last round", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run()
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  out <- dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1)
  expect_equal(seen$files, c("t1_a.log", "t1_a.trees", "t1_a.xml.state",
    "t1_b.log", "t1_b.trees", "t1_b.xml.state"))
  # the growing files are cut to the size at the checkpoint (2 lines)
  expect_equal(seen$content[["t1_a.log"]], c("round", "round"))
  # a tree log gets back its closing line
  expect_equal(seen$content[["t1_a.trees"]], c("round", "round", "End;"))
  # the small files are the copies of the checkpoint
  expect_equal(seen$content[["t1_a.xml.state"]], c("round", "round"))
  expect_equal(seen$content[["t1_b.log"]], c("round", "round"))
})

test_that("the files in dir are not changed, and the temporary folder is deleted", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run()
  files <- function() {
    f <- list.files(s$dir, recursive=TRUE, full.names=TRUE)
    stats::setNames(unname(tools::md5sum(f)), f)
  }
  before <- files()
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1)
  expect_identical(files(), before)
  expect_equal(length(readLines(file.path(s$dir, "t1_a.log"))), 3)
  expect_false(identical(seen$dir, s$dir))
  expect_false(dir.exists(seen$dir))
})

test_that("the temporary folder is also deleted if readBEAST fails", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run()
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen, fail=TRUE), .package="dowser")
  expect_error(dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1),
    "readBEAST failed")
  expect_false(dir.exists(seen$dir))
})

test_that("a tree log that ended with End; gets only one", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- iterate_setup("a", function(cid, n) 1000)
  s$work(s$clones[1, ], NULL)
  result <- s$work(s$clones[1, ], "a")
  writeBin(charToRaw("tree 1;\nEnd;"), file.path(s$dir, "t1_a.trees"))
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 2L, result)
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1)
  expect_equal(seen$content[["t1_a.trees"]], c("tree 1;", "End;"))
})

test_that("the tree log with traits is rebuilt in the same way", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- iterate_setup("a", function(cid, n) 1000)
  s$work(s$clones[1, ], NULL)
  result <- s$work(s$clones[1, ], "a")
  writeBin(charToRaw("tree 1;\nEnd;"), file.path(s$dir, "t1_a_tree_with_trait.trees"))
  dowser:::writeIterateCheckpoint(s$dir, "t1", "a", "a", 2L, result)
  cat("tree 2;\n", file=file.path(s$dir, "t1_a_tree_with_trait.trees"), append=TRUE)
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", trait="location", quiet=1)
  expect_equal(seen$content[["t1_a_tree_with_trait.trees"]], c("tree 1;", "End;"))
})

test_that("the arguments are passed to readBEAST, and the default posterior is none", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run("a")
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1)
  expect_equal(seen$args$posterior, "none")
  expect_equal(seen$args$id, "t1")
  expect_equal(seen$args$beast, "beast_bin")
  expect_equal(seen$args$burnin, 10)
  dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "my_beast", burnin=20, trait="location",
    nproc=2, quiet=1, posterior="trees", asr=TRUE, low_ram=FALSE, trim_ids=TRUE)
  expect_equal(seen$args, list(id="t1", beast="my_beast", burnin=20, trait="location", nproc=2,
    quiet=1, posterior="trees", asr=TRUE, low_ram=FALSE, trim_ids=TRUE, clone_ids="a"))
})

test_that("the result has below_ESS and the round of each checkpoint", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run()
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  out <- dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1)
  expect_equal(out$clone_id, c("a", "b"))
  expect_equal(out$below_ESS, c(0L, 0L))
  expect_equal(out$checkpoint_round, c(2L, 2L))
  out <- dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1, ess_cutoff=5000)
  expect_equal(out$below_ESS, c(2L, 2L))
})

test_that("clones with different rounds keep their own round", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- iterate_setup(c("a", "b"), function(cid, n) 1000)
  for(cid in c("a", "b")){
    i <- match(cid, c("a", "b"))
    result <- s$work(s$clones[i, ], NULL)
    n_rounds <- if(cid == "a") 1L else 3L
    for(r in seq_len(n_rounds - 1)) result <- s$work(s$clones[i, ], cid)
    dowser:::writeIterateCheckpoint(s$dir, "t1", c("a", "b"), cid, n_rounds, result)
  }
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  out <- dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1)
  expect_equal(out$checkpoint_round, c(1L, 3L))
  expect_equal(seen$content[["t1_a.log"]], "round")
  expect_equal(seen$content[["t1_b.log"]], c("round", "round", "round"))
})

test_that("a clone without a good checkpoint is left out, with a warning", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run()
  unlink(file.path(s$dir, "checkpoint_t1", "clone_b"), recursive=TRUE)
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  expect_warning(out <- dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin",
    quiet=1), "clone\\(s\\) b")
  expect_equal(out$clone_id, "a")
  expect_equal(seen$args$clone_ids, "a")
  expect_equal(seen$files, c("t1_a.log", "t1_a.trees", "t1_a.xml.state"))
})

test_that("a clone whose growing file is shorter than at the checkpoint is left out", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run()
  writeLines("x", file.path(s$dir, "t1_a.log"))
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  expect_warning(out <- dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin",
    quiet=1), "clone\\(s\\) a")
  expect_equal(out$clone_id, "b")
})

test_that("it stops if no clone has a good checkpoint", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- iterate_setup(c("a", "b"), function(cid, n) 1000)
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  expect_error(suppressWarnings(dowser::readFromLastCheckpoint(s$clones, s$dir, "t1",
    "beast_bin", quiet=1)), "No clone has a good checkpoint")
  expect_null(seen$dir)
})

test_that("it says the round of each clone unless quiet", {
  skip_if_not_installed("testthat", "3.2.0")
  s <- unfinished_run()
  seen <- new.env()
  testthat::local_mocked_bindings(readBEAST=spy_readBEAST(seen), .package="dowser")
  txt <- capture.output(dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin"))
  expect_true(any(grepl("Clone a : results as of round 2", txt, fixed=TRUE)))
  expect_true(any(grepl("Clone b : results as of round 2", txt, fixed=TRUE)))
  txt <- capture.output(dowser::readFromLastCheckpoint(s$clones, s$dir, "t1", "beast_bin", quiet=1))
  expect_false(any(grepl("results as of round", txt, fixed=TRUE)))
})
