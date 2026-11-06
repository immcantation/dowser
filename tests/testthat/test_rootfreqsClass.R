# Tests for the choice of the sequence likelihood class that has root frequencies
# (tyche.evolution.likelihood.RootFreqsTreeLikelihood for TyCHE >= 1.0.1, otherwise
# rootfreqs.TreeLikelihood).
#
# Part 1: the pure functions (extract_tyche_version, parse_tyche_version, choose_rootfreqs_class).
# Part 2: resolve_rootfreqs_class() and buildBeast(), with fake "beast" executables (shell
#         scripts that print a chosen `-version` output). No real BEAST is used. The fake
#         executables need a shell, so Part 2 is skipped on Windows.
#
# These functions are not exported, so they are called directly by name here.

# nolint start

TYCHE_CLASS  <- "tyche.evolution.likelihood.RootFreqsTreeLikelihood"
LEGACY_CLASS <- "rootfreqs.TreeLikelihood"

# ---- Part 1: pure functions --------------------------------------------------

beast27_output <- c(
  "BEAST v2.7.7", "---", "BEAST.app v2.7.7", "BEAST.base v2.7.7",
  "TyCHE v-cb8e6a7-built-10221605", "rootfreqs v0.0.2", "BEAST_CLASSIC v1.6.3",
  "---", "Java version 17.0.3")

beast3_output <- c(
  "Picked up JAVA_TOOL_OPTIONS: -Xmx2g",
  "BEAST modules in boot layer: [beast.base, beast.pkgmgmt]",
  "Loading package BEAST.base v2.8.0 from /home/u/.beast/2.8/BEAST.base",
  "Loading package TyCHE v9.9.9 from /home/u/.beast/2.8/TyCHE",
  "BEAST v2.8.0", "---", "BEAST.base v2.8.0", "TyCHE v1.0.1-ab12cd3-dirty-built-10271545",
  "BEAST_CLASSIC v1.7.0", "---", "Java version 25")

test_that("extract_tyche_version(): takes the rest of the line that starts with 'TyCHE v'", {
  expect_identical(extract_tyche_version(beast27_output), "-cb8e6a7-built-10221605")
  # BEAST 3: only the line that starts with 'TyCHE v', not the 'Loading package TyCHE v9.9.9' line
  expect_identical(extract_tyche_version(beast3_output), "1.0.1-ab12cd3-dirty-built-10271545")
})

test_that("extract_tyche_version(): NULL if there is no TyCHE line; copes with CR and spaces", {
  expect_null(extract_tyche_version(c("BEAST v2.7.7", "---", "BEAST.base v2.7.7", "---")))
  expect_null(extract_tyche_version(character(0)))
  expect_null(extract_tyche_version("Loading package TyCHE v1.2.3 from /x"))  # substring only
  expect_null(extract_tyche_version("TyCHE 1.0.1"))                           # no 'v'
  expect_identical(extract_tyche_version(c("BEAST v2.7.7", "TyCHE v1.0.1\r")), "1.0.1")
  expect_identical(extract_tyche_version("  TyCHE v1.0.1  "), "1.0.1")
  # the first TyCHE line is used
  expect_identical(extract_tyche_version(c("TyCHE v1.0.1", "TyCHE v0.0.1")), "1.0.1")
})

test_that("parse_tyche_version(): leading MAJOR.MINOR.PATCH, rest ignored, missing parts are 0", {
  expect_equal(parse_tyche_version("1.0.1"), c(1, 0, 1))
  expect_equal(parse_tyche_version("1.0.1-ab12cd3-built-10271545"), c(1, 0, 1))
  expect_equal(parse_tyche_version("1.0.1-ab12cd3-dirty-built-10271545"), c(1, 0, 1))
  expect_equal(parse_tyche_version("1.0.10"), c(1, 0, 10))
  expect_equal(parse_tyche_version("1.1"), c(1, 1, 0))
  expect_equal(parse_tyche_version("2"), c(2, 0, 0))
  expect_equal(parse_tyche_version("0.0.11"), c(0, 0, 11))
  expect_equal(parse_tyche_version("1.0.1.7"), c(1, 0, 1))
  # no leading number
  expect_null(parse_tyche_version("-cb8e6a7-built-10221605"))
  expect_null(parse_tyche_version(""))
  expect_null(parse_tyche_version("abc"))
  expect_null(parse_tyche_version(NA_character_))
  expect_null(parse_tyche_version(NULL))
})

test_that("choose_rootfreqs_class(): the table of installed versions", {
  to_tyche <- c("1.0.1", "1.0.1-ab12cd3-built-10271545", "1.0.1-ab12cd3-dirty-built-10271545",
    "1.0.2", "1.0.2-ab12cd3-built-10271545", "1.0.10", "1.1.0", "2.0.0", "1.1", "10.0.0")
  to_legacy <- c("1.0.0", "1.0.0-ab12cd3-built-10271545", "0.0.11", "0.9.9", "0.99.99", "0")
  for (v in to_tyche) expect_identical(choose_rootfreqs_class(v), TYCHE_CLASS, info = v)
  for (v in to_legacy) expect_identical(choose_rootfreqs_class(v), LEGACY_CLASS, info = v)
  # cannot decide: legacy class
  expect_identical(choose_rootfreqs_class("-cb8e6a7-built-10221605"), LEGACY_CLASS)
  expect_identical(choose_rootfreqs_class(NULL), LEGACY_CLASS)
  expect_identical(choose_rootfreqs_class(NA_character_), LEGACY_CLASS)
})

test_that("choose_rootfreqs_class(): numbers are compared as numbers, not as text", {
  # as text "1.0.10" < "1.0.2" and "10.0.0" < "2.0.0"
  expect_identical(choose_rootfreqs_class("1.0.10"), TYCHE_CLASS)
  expect_identical(choose_rootfreqs_class("0.10.0"), LEGACY_CLASS)  # 0.10.0 is below 1.0.1
  expect_identical(choose_rootfreqs_class("0.0.100"), LEGACY_CLASS)
  expect_identical(choose_rootfreqs_class("1.0.100"), TYCHE_CLASS)
})

test_that("check_rootfreqs_class(): accepts class names, stops for anything else", {
  expect_identical(check_rootfreqs_class(TYCHE_CLASS), TYCHE_CLASS)
  expect_identical(check_rootfreqs_class(LEGACY_CLASS), LEGACY_CLASS)
  expect_error(check_rootfreqs_class(""), "single BEAST class name")
  expect_error(check_rootfreqs_class(NULL), "single BEAST class name")
  expect_error(check_rootfreqs_class(c("a.B", "c.D")), "single BEAST class name")
  expect_error(check_rootfreqs_class('a"b'), "single BEAST class name")
})

# ---- Part 2: fake BEAST executables ------------------------------------------

# A directory with fake beast, treeannotator and loganalyser. `beast -version` prints
# `version_lines` (and a JAVA_TOOL_OPTIONS line on stderr) and exits with `status`.
# Each run of fake beast adds a line to the file "calls" in the directory.
.make_fake_beast <- function(version_lines, status = 0) {
  dir <- tempfile("fakebeast_")
  dir.create(dir)
  body <- c("#!/bin/sh",
    'echo called >> "$(dirname "$0")/calls"',
    if (length(version_lines) > 0) paste0("echo '", version_lines, "'"),
    'echo "Picked up JAVA_TOOL_OPTIONS: -Xmx1g" >&2',
    paste0("exit ", status))
  for (f in c("beast", "treeannotator", "loganalyser")) {
    writeLines(body, file.path(dir, f))
    Sys.chmod(file.path(dir, f), "755")
  }
  dir
}

.n_calls <- function(dir) {
  f <- file.path(dir, "calls")
  if (!file.exists(f)) return(0L)
  length(readLines(f))
}

.version_output <- function(tyche_line) {
  c("BEAST v2.7.7", "---", "BEAST.app v2.7.7", "BEAST.base v2.7.7", tyche_line,
    "rootfreqs v0.0.2", "BEAST_CLASSIC v1.6.3", "---", "Java version 17.0.3")
}

tmpl_rf_path <- tempfile(fileext = ".xml")
writeLines(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '<distribution id="treeLikelihood.${CLONE}" spec="TreeLikelihood" data="@${CLONE}" tree="@Tree.t:${CLONE}">',
  '${ROOTFREQS}',
  '</distribution>',
  '</run>',
  '</beast>'), tmpl_rf_path)

tmpl_no_rf_path <- tempfile(fileext = ".xml")
writeLines(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '<distribution id="treeLikelihood.${CLONE}" spec="TreeLikelihood" data="@${CLONE}" tree="@Tree.t:${CLONE}">',
  '</distribution>',
  '</run>',
  '</beast>'), tmpl_no_rf_path)

.rf_clone <- function(clone_id = "1") {
  methods::new("airrClone",
    data = data.frame(
      sequence_id = c("seq1", "seq2", "seq3"),
      sequence = c("AAAAAAAA", "CCCCCCCC", "GGGGTTTT"),
      collection_date = c(2000, 2005, 2010),
      celltype = c("GC", "GC", "MBC"),
      stringsAsFactors = FALSE),
    clone = clone_id, germline = "AACCGGTT", lgermline = "", hlgermline = "",
    v_gene = "IGHV1", j_gene = "IGHJ1", junc_len = 8, locus = "IGH",
    region = character(0), phylo_seq = "sequence", numbers = numeric(0))
}

# buildBeast() stops after the XML is written, and this returns the paths of the XML files
.xml_paths_from_buildBeast <- function(beast, template, clones = list(.rf_clone("1")), ...) {
  out <- tempfile("rfout_")
  dir.create(out)
  stop_after_xml <- function(paths) {
    stop(structure(class = c("rf_xml_only", "error", "condition"),
      list(message = "XML written", call = NULL, paths = paths)))
  }
  tryCatch(
    buildBeast(clones, beast = beast, time = "collection_date", template = template,
      dir = out, id = "run", trait = "celltype", quiet = 1, xml_check = stop_after_xml, ...),
    rf_xml_only = function(cnd) cnd$paths)
}

test_that("resolve_rootfreqs_class(): picks the class from `beast -version`", {
  skip_on_os("windows")
  cases <- list(
    list("TyCHE v1.0.1", TYCHE_CLASS),
    list("TyCHE v1.0.1-ab12cd3-dirty-built-10271545", TYCHE_CLASS),
    list("TyCHE v1.0.10", TYCHE_CLASS),
    list("TyCHE v1.0.0-ab12cd3-built-10271545", LEGACY_CLASS),
    list("TyCHE v0.0.11", LEGACY_CLASS))
  for (case in cases) {
    beast <- .make_fake_beast(.version_output(case[[1]]))
    got <- suppressMessages(resolve_rootfreqs_class(beast, tmpl_rf_path, TRUE, NULL, 0))
    expect_identical(got, case[[2]], info = case[[1]])
    expect_equal(.n_calls(beast), 1)
  }
})

test_that("resolve_rootfreqs_class(): cannot decide -> legacy class and a message that says why", {
  skip_on_os("windows")
  # an old dev build: no leading number
  beast <- .make_fake_beast(.version_output("TyCHE v-cb8e6a7-built-10221605"))
  expect_message(
    got <- resolve_rootfreqs_class(beast, tmpl_rf_path, TRUE, NULL, 1),
    "Could not decide.*does not start with a number.*rootfreqs.TreeLikelihood")
  expect_identical(got, LEGACY_CLASS)

  # no TyCHE line at all
  beast <- .make_fake_beast(c("BEAST v2.7.7", "---", "BEAST.base v2.7.7", "---"))
  expect_message(
    got <- resolve_rootfreqs_class(beast, tmpl_rf_path, TRUE, NULL, 1),
    "Could not decide.*no line starting with 'TyCHE v'")
  expect_identical(got, LEGACY_CLASS)

  # BEAST that prints nothing and fails
  beast <- .make_fake_beast(character(0), status = 1)
  expect_message(
    got <- resolve_rootfreqs_class(beast, tmpl_rf_path, TRUE, NULL, 1),
    "Could not decide")
  expect_identical(got, LEGACY_CLASS)

  # no BEAST at that path
  expect_message(
    got <- resolve_rootfreqs_class(file.path(tempdir(), "no_such_beast_dir"), tmpl_rf_path, TRUE, NULL, 1),
    "Could not decide")
  expect_identical(got, LEGACY_CLASS)
})

test_that("resolve_rootfreqs_class(): BEAST is not run if the class is given or not needed", {
  skip_on_os("windows")
  beast <- .make_fake_beast(.version_output("TyCHE v1.0.1"))
  # the caller gives the class (even a class the version would not pick)
  expect_identical(resolve_rootfreqs_class(beast, tmpl_rf_path, TRUE, LEGACY_CLASS, 0), LEGACY_CLASS)
  expect_identical(resolve_rootfreqs_class(beast, tmpl_rf_path, TRUE, TYCHE_CLASS, 0), TYCHE_CLASS)
  expect_error(resolve_rootfreqs_class(beast, tmpl_rf_path, TRUE, "bad name", 0), "single BEAST class name")
  # germline is not the root
  expect_null(resolve_rootfreqs_class(beast, tmpl_rf_path, FALSE, NULL, 0))
  # template without ${ROOTFREQS}
  expect_null(resolve_rootfreqs_class(beast, tmpl_no_rf_path, TRUE, NULL, 0))
  # template file that does not exist (write_clone_to_xml reports it), and no beast
  expect_null(resolve_rootfreqs_class(beast, file.path(tempdir(), "no_template.xml"), TRUE, NULL, 0))
  expect_null(resolve_rootfreqs_class(NULL, tmpl_rf_path, TRUE, NULL, 0))
  expect_equal(.n_calls(beast), 0)
})

test_that("buildBeast(): writes the class of the installed TyCHE into the XML, BEAST runs once", {
  skip_on_os("windows")
  # new TyCHE
  beast <- .make_fake_beast(.version_output("TyCHE v1.0.2-ab12cd3-built-10271545"))
  paths <- .xml_paths_from_buildBeast(beast, tmpl_rf_path)
  content <- paste(readLines(paths[1]), collapse = "\n")
  expect_true(grepl(paste0('spec="', TYCHE_CLASS, '"'), content, fixed = TRUE))
  expect_false(grepl("rootfreqs.TreeLikelihood", content, fixed = TRUE))
  expect_true(grepl('<rootfreqseq id="seq_Root" spec="Sequence" taxon="Root"', content, fixed = TRUE))
  expect_equal(.n_calls(beast), 1)

  # old TyCHE: today's behavior
  beast <- .make_fake_beast(.version_output("TyCHE v0.0.11"))
  paths <- .xml_paths_from_buildBeast(beast, tmpl_rf_path)
  content <- paste(readLines(paths[1]), collapse = "\n")
  expect_true(grepl('spec="rootfreqs.TreeLikelihood"', content, fixed = TRUE))
  expect_false(grepl("RootFreqsTreeLikelihood", content, fixed = TRUE))
  expect_equal(.n_calls(beast), 1)
})

test_that("buildBeast(): version is read once per call, not once per clone", {
  skip_on_os("windows")
  beast <- .make_fake_beast(.version_output("TyCHE v1.0.1"))
  clones <- list(.rf_clone("1"), .rf_clone("2"), .rf_clone("3"))
  paths <- .xml_paths_from_buildBeast(beast, tmpl_rf_path, clones = clones)
  expect_length(paths, 3)
  for (p in paths) {
    expect_true(grepl(TYCHE_CLASS, paste(readLines(p), collapse = "\n"), fixed = TRUE))
  }
  expect_equal(.n_calls(beast), 1)
})

test_that("buildBeast(): rootfreqs_class given -> BEAST is not asked for its version", {
  skip_on_os("windows")
  beast <- .make_fake_beast(.version_output("TyCHE v1.0.1"))   # would pick the tyche class
  paths <- .xml_paths_from_buildBeast(beast, tmpl_rf_path, rootfreqs_class = LEGACY_CLASS)
  content <- paste(readLines(paths[1]), collapse = "\n")
  expect_true(grepl('spec="rootfreqs.TreeLikelihood"', content, fixed = TRUE))
  expect_equal(.n_calls(beast), 0)
})

test_that("buildBeast(): BEAST is not asked for its version if the XML does not need the class", {
  skip_on_os("windows")
  beast <- .make_fake_beast(.version_output("TyCHE v1.0.1"))
  # template has no ${ROOTFREQS}
  paths <- .xml_paths_from_buildBeast(beast, tmpl_no_rf_path)
  content <- paste(readLines(paths[1]), collapse = "\n")
  expect_true(grepl('spec="TreeLikelihood"', content, fixed = TRUE))
  # germline not used as root
  paths <- .xml_paths_from_buildBeast(beast, tmpl_rf_path, include_germline = FALSE)
  content <- paste(readLines(paths[1]), collapse = "\n")
  expect_true(grepl('spec="TreeLikelihood"', content, fixed = TRUE))
  expect_false(grepl("rootfreqseq", content, fixed = TRUE))
  expect_equal(.n_calls(beast), 0)
})

# nolint end
