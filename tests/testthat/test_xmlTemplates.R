# Tests for the code that writes the BEAST XML files, and for the templates that it uses.
#
# Part 1: write_clone_to_xml() / write_clones_to_xmls(), and the pure XML-builder helpers that
# they call (create_alignment, create_traitset, create_root_freqs, create_MRCA_prior_observed and
# _germline, create_height_prior, create_max_height_prior, create_starting_tree).
# Part 2 (at the end): the dowser-template-info convention: the comment-block parser
# (.parseTemplateInfo), the built-in structural checklist (.structuralChecklist), the exported
# describeTemplate() with its full dump and its single-variable (var=) lookup, and the error
# hook (.templateVarLookup) that write_clone_to_xml() uses.
#
# Most of these functions are not exported, so they are called directly by name here.
# No BEAST2 binary is invoked anywhere in this file: everything is text manipulation on toy
# templates and a hand-built toy airrClone object.
#
# One source-level quirk is locked in as documented behavior, with a comment at its test: the
# date of the germline is never added to the tip dates.

# ---- Test utilities --------------------------------------------------------

# nolint start
.write_template <- function(lines) {
  path <- tempfile(fileext = ".xml")
  writeLines(lines, path)
  path
}

.read_all <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

.extract <- function(text, pattern) {
  m <- regmatches(text, regexpr(pattern, text, perl = TRUE))
  if (length(m) == 0 || identical(m, "")) return(NA_character_)
  m
}

.count_matches <- function(text, pattern) {
  length(regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1]])
}

# mirrors the codeMap-building snippet inside write_clones_to_xmls(), so
# expected values are derived from the same primitive paste()/paste0()
# calls rather than hand-typed literals
.build_trait_data_type <- function(trait_list) {
  codeMap <- paste(trait_list, 0:(length(trait_list) - 1), collapse = ",\n", sep = "=")
  codeMap <- paste0(codeMap, ",\n? = ", paste0(0:(length(trait_list) - 1), collapse = " "))
  paste0('<userDataType id="traitDataType.newTrait" spec="beast.base.evolution.datatype.UserDataType" codeMap="',
    codeMap, '" codelength="-1" states="', length(trait_list), '"/>')
}

# ---- Toy airrClone fixture --------------------------------------------------
#
# 3 short sequences (8nt, no ambiguous bases so germline/root-freq behavior
# is predictable), a numeric date column, and a 2-state categorical trait
# column. Built directly via new("airrClone", ...) rather than via
# makeAirrClone()/formatClones() so every field is fully known and
# hand-verifiable -- this test file is only exercising the XML string
# substitution logic, not clone construction.
.new_toy_clone <- function(clone_id = "1", germline = "AACCGGTT") {
  methods::new("airrClone",
    data = data.frame(
      sequence_id = c("seq1", "seq2", "seq3"),
      sequence = c("AAAAAAAA", "CCCCCCCC", "GGGGTTTT"),
      collection_date = c(2000, 2005, 2010),
      celltype = c("GC", "GC", "MBC"),
      stringsAsFactors = FALSE
    ),
    clone = clone_id,
    germline = germline,
    lgermline = "",
    hlgermline = "",
    v_gene = "IGHV1",
    j_gene = "IGHJ1",
    junc_len = 8,
    locus = "IGH",
    region = character(0),
    phylo_seq = "sequence",
    numbers = numeric(0)
  )
}

clone1 <- .new_toy_clone(clone_id = "1", germline = "AACCGGTT")
trait_list_toy <- c("GC", "MBC") # sort(unique(clone1@data$celltype))
trait_data_type_toy <- .build_trait_data_type(trait_list_toy)

# write_clone_to_xml() prints the nucleotide frequencies when a template has ${EMP_EQFREQS},
# and the trait list before the "unknown state" error. A run with nproc=1 relies on the first
# one to show the frequencies, so the code keeps it, and the tests of this file hide it.
# write_clone_to_xml and write_clones_to_xmls here are the functions of the package with their
# printed output removed (warnings and errors are not changed, so expect_warning() and
# expect_error() work as before). A call from inside the package is not changed.
.without_output <- function(f){
  force(f)    # f must be the function of the package, not the new one that has the same name
  function(...){
    out <- NULL
    utils::capture.output(out <- f(...))
    out
  }
}
write_clone_to_xml <- .without_output(write_clone_to_xml)
write_clones_to_xmls <- .without_output(write_clones_to_xmls)

# ---- Toy templates ----------------------------------------------------------
#
# A: GermlineRootTree-style tree, no optional placeholders.
# B: GermlineRootTree-style tree, every optional placeholder present.
# C: same optional placeholders as B, but no GermlineRootTree marker text.
# D: generic catch-all free-form placeholders + <init> (for replacements()).
# TC: minimal template with a nodeTypes parameter (for tree_states).
# WRAP: minimal template used only for the write_clones_to_xmls() wrapper
#       tests (DATA/CLONE/TRAIT/TRAIT_NAME/init only, so multi-clone
#       behavior isn't obscured by unrelated warnings from other
#       placeholders).
# Broken: E1 has no <init> tag at all.

tmpl_grt_minimal_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<state>',
  '<tree id="Tree.t:${CLONE}" spec="tyche.evolution.tree.GermlineRootTree" name="stateNode">',
  '</tree>',
  '</state>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '${OPERATORS}',
  '<logger id="tracelog" fileName="${CLONE}.log" logEvery="1000">',
  '</logger>',
  '<logger id="treeWithTraitLogger" fileName="${CLONE}.trees" logEvery="1000" mode="tree">',
  '</logger>',
  '</run>',
  '</beast>'
))

.grt_full_lines <- function(tree_spec) c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<parameter id="nodeTypes" dimension="${NODES}">0</parameter>',
  '<state>',
  paste0('<tree id="Tree.t:${CLONE}" spec="', tree_spec, '" name="stateNode">'),
  '${DATE}',
  '<taxonset idref="TaxonSet.${CLONE}"/>',
  '</tree>',
  '</state>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '<distribution id="prior" spec="CompoundDistribution">',
  '${MRCA}',
  '<distribution id="typedTreeLikelihood.type" spec="tyche.evolution.likelihood.AncestralTypeLikelihood" tag="${TRAIT_NAME}" tree="@Tree.t:${CLONE}">',
  '${TRAIT}',
  '${ROOT_TYPE_PROBABILITIES}',
  '</distribution>',
  '</distribution>',
  '<distribution id="treeLikelihood.${CLONE}" spec="TreeLikelihood" data="@${CLONE}" tree="@Tree.t:${CLONE}">',
  '${ROOTFREQS}',
  '<parameter id="freqParameter.s:${CLONE}" dimension="4">${EMP_EQFREQS}</parameter>',
  '</distribution>',
  '<parameter id="nodeTypesInit" value="${NODES_TYPE_INIT}"/>',
  '${OPERATORS}',
  '<logger id="tracelog" fileName="${CLONE}.log" logEvery="1000">',
  '</logger>',
  '<logger id="treeWithTraitLogger" fileName="${CLONE}.trees" logEvery="1000" mode="tree">',
  '</logger>',
  '</run>',
  '</beast>'
)

tmpl_grt_full_path    <- .write_template(.grt_full_lines("tyche.evolution.tree.GermlineRootTree"))
tmpl_no_grt_full_path <- .write_template(.grt_full_lines("beast.base.evolution.tree.Tree"))

tmpl_catchall_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<parameter id="fooParam.${CLONE}" value="${FOO}"/>',
  '<parameter id="barParam.${CLONE}" value="${BAR}"/>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '</run>',
  '</beast>'
))

tmpl_grt_minimal_tc_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<parameter id=\'nodeTypes\' spec=\'parameter.IntegerParameter\' dimension="${NODES}" value="0"/>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '</run>',
  '</beast>'
))

tmpl_grt_minimal_no_tc_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '</run>',
  '</beast>'
))

tmpl_old_trait_style_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<distribution id="typedTreeLikelihood.type" spec="tyche.evolution.likelihood.AncestralTypeLikelihood" tag="location" tree="@Tree.t:${CLONE}">',
  '${TRAIT}',
  '</distribution>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '</run>',
  '</beast>'
))

tmpl_broken_no_init_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<state></state>',
  '</run>',
  '</beast>'
))

tmpl_wrap_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<tree id="Tree.t:${CLONE}" spec="beast.base.evolution.tree.Tree" name="stateNode">',
  '<taxonset idref="TaxonSet.${CLONE}"/>',
  '</tree>',
  '<distribution id="typedTreeLikelihood.type" spec="tyche.evolution.likelihood.AncestralTypeLikelihood" tag="${TRAIT_NAME}" tree="@Tree.t:${CLONE}">',
  '${TRAIT}',
  '</distribution>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '</run>',
  '</beast>'
))

# ---- ${DATA} / ${CLONE} / chainLength / logEvery ---------------------------

test_that("core substitution: ${DATA}, ${CLONE}, chainLength, logEvery", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
    log_every_trait = 10, include_germline_as_tip = FALSE)

  expect_equal(path, paste0(outfile, "_1.xml"))
  expect_true(file.exists(path))
  content <- .read_all(path)

  # ${CLONE} substituted everywhere as id_clone
  expect_true(grepl('Tree.t:run_1', content, fixed = TRUE))
  expect_true(grepl('taxa="@run_1"', content, fixed = TRUE))

  # ${DATA}: one <sequence> per tip, no Germline tip since include_germline_as_tip=FALSE
  expect_true(grepl('<sequence id="seq_seq1" spec="Sequence" taxon="seq1" totalcount="4" value="AAAAAAAA" />', content, fixed = TRUE))
  expect_true(grepl('<sequence id="seq_seq2" spec="Sequence" taxon="seq2" totalcount="4" value="CCCCCCCC" />', content, fixed = TRUE))
  expect_true(grepl('<sequence id="seq_seq3" spec="Sequence" taxon="seq3" totalcount="4" value="GGGGTTTT" />', content, fixed = TRUE))

  expect_false(grepl('taxon="Germline"', content, fixed = TRUE))
  expect_false(grepl('id="seq_Germline"', content, fixed = TRUE))

  expect_true(grepl('<data id="run_1" spec="Alignment" name="alignment">', content, fixed = TRUE))
  expect_true(grepl('<taxa id="TaxonSet.run_1" spec="TaxonSet">', content, fixed = TRUE))

  # chainLength overwritten with mcmc_length
  expect_true(grepl('chainLength="5000"', content, fixed = TRUE))

  # logEvery: treeWithTraitLogger gets log_every * log_every_trait, others get log_every
  expect_true(grepl('logEvery="1000"', content, fixed = TRUE)) # 100 * 10, trait logger
  expect_true(grepl('logEvery="100"', content, fixed = TRUE))  # tracelog

  # no unfilled placeholders remain
  expect_false(grepl("\\$\\{", content))
})

test_that("core substitution: ${DATA} includes Germline tip when requested", {
  outfile <- tempfile()
# tmpl_grt_minimal_path contains the GermlineRootTree marker, so
  # include_germline_as_tip=TRUE also triggers the OPERATORS skip-warning
  # tested separately below; expected here too, not the focus of this test.
  expect_warning(
    path <- write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = TRUE),
    "GermlineRootTree"
  )
  content <- .read_all(path)

  expect_true(grepl('<sequence id="seq_Germline" spec="Sequence" taxon="Germline" totalcount="4" value="AACCGGTT" />', content, fixed = TRUE))
})

# ---- ${OPERATORS}: GermlineRootTree skip-with-warning vs. add operator ----

test_that("${OPERATORS}: GermlineRootTree marker skips TipDatesRandomWalker with a warning", {
  outfile <- tempfile()
  expect_warning(
    path <- write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = TRUE),
    "GermlineRootTree"
  )
  content <- .read_all(path)
  expect_false(grepl("TipDatesRandomWalker", content, fixed = TRUE))
})

test_that("${OPERATORS}: without GermlineRootTree marker, adds TipDatesRandomWalker with no warning", {
  outfile <- tempfile()
  path <- withCallingHandlers(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_no_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
      include_germline_as_tip = TRUE, include_germline_as_root = FALSE),
    warning = function(w) fail(paste("unexpected warning:", conditionMessage(w)))
  )
  content <- .read_all(path)
  expect_true(grepl(
    '<operator id="TipDatesRandomWalker.01" windowSize="1" spec="beast.base.evolution.operator.TipDatesRandomWalker" taxonset="@germSet" tree="@Tree.t:run_1" weight="1.0"/>',
    content, fixed = TRUE))
})

test_that("${OPERATORS}: empty (no warning, no operator) when include_germline_as_tip is FALSE", {
  outfile <- tempfile()
  path <- withCallingHandlers(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = FALSE),
    warning = function(w) fail(paste("unexpected warning:", conditionMessage(w)))
  )
  content <- .read_all(path)
  expect_false(grepl("TipDatesRandomWalker", content, fixed = TRUE))
})

# ---- ${DATE} ----------------------------------------------------------------

test_that("${DATE}: date=NULL raises an error when ${DATE} is present", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100),
    "Date argument is NULL"
  )
})

test_that("${DATE}: traitSet is built from the date column, and NEVER includes a Germline row", {
  outfile <- tempfile()
  # tmpl_grt_full_path contains the GermlineRootTree marker; see note above
  expect_warning(
    path <- write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
      include_germline_as_tip = TRUE, root_trait = NULL),
    "GermlineRootTree"
  )
  content <- .read_all(path)

  date_block <- .extract(content, '<trait id="date:run_1"[\\s\\S]*?</trait>')
  expect_false(is.na(date_block))
  expect_true(grepl('spec="beast.base.evolution.tree.TraitSet" taxa="@TaxonSet.run_1" traitname="date"', date_block, fixed = TRUE))
  expect_true(grepl("seq1=2000", date_block, fixed = TRUE))
  expect_true(grepl("seq2=2005", date_block, fixed = TRUE))
  expect_true(grepl("seq3=2010", date_block, fixed = TRUE))
  # documented quirk: germline date is never added, even though
  # include_germline_as_tip = TRUE here and the germline IS a tip elsewhere
  expect_false(grepl("Germline", date_block, fixed = TRUE))
})

# ---- ${TRAIT} / ${TRAIT_NAME} -----------------------------------------------

test_that("${TRAIT}: trait=NULL raises an error when ${TRAIT} is present", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date"),
    "Trait argument is NULL"
  )
})

test_that("${TRAIT}: trait given without trait_data_type raises an error", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date", trait = "celltype"),
    "no trait data type"
  )
})

test_that("${TRAIT}/${TRAIT_NAME}: traitSet built from trait column, DOES include a Germline row", {
  outfile <- tempfile()
  # tmpl_grt_full_path contains the GermlineRootTree marker; see note above
  expect_warning(
    path <- write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
      include_germline_as_tip = TRUE),
    "GermlineRootTree"
  )
  content <- .read_all(path)

  # ${TRAIT_NAME} filled with the raw column name
  expect_true(grepl('tag="celltype"', content, fixed = TRUE))

  trait_block <- .extract(content, '<traitSet id="newTrait:run_1"[\\s\\S]*?</traitSet>')
  expect_false(is.na(trait_block))
  expect_true(grepl('spec="beast.base.evolution.tree.TraitSet" taxa="@TaxonSet.run_1" traitname="newTrait"', trait_block, fixed = TRUE))
  expect_true(grepl("seq1=GC", trait_block, fixed = TRUE))
  expect_true(grepl("seq2=GC", trait_block, fixed = TRUE))
  expect_true(grepl("seq3=MBC", trait_block, fixed = TRUE))
  # unlike ${DATE}, ${TRAIT} DOES include the germline (default value '?')
  expect_true(grepl("Germline=\\?", trait_block))

  # trait_data_type appended right after the traitSet
  expect_true(grepl(trait_data_type_toy, content, fixed = TRUE))
})

test_that("old-style templates (no ${TRAIT_NAME}) fall back to rewriting tag=\"location\"", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_old_trait_style_path, mcmc_length = 5000, log_every = 100,
    trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE)
  content <- .read_all(path)

  expect_true(grepl('tag="celltype"', content, fixed = TRUE))
  expect_false(grepl('tag="location"', content, fixed = TRUE))
})

# ---- ${NODES} ---------------------------------------------------------------

test_that("${NODES}: node count is 2*(nrow(data)+1)-1, germline counted as a tip", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE)
  content <- .read_all(path)
  # nrow(clone1@data) == 3 -> tips = 4 -> nodes = 7
  expect_true(grepl('dimension="7"', content, fixed = TRUE))
})

# ---- ${MRCA} -----------------------------------------------------------------

test_that("${MRCA}: empty when include_germline_as_tip=FALSE and no start/max dates", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE)
  content <- .read_all(path)
  # the enclosing <distribution id="prior"> should have nothing between
  # it and the next known tag when MRCA is empty
  expect_false(grepl("obs.prior", content, fixed = TRUE))
  expect_false(grepl("germ1.prior", content, fixed = TRUE))
  expect_false(grepl("height.prior", content, fixed = TRUE))
})

test_that("${MRCA}: observed + germline MRCA priors present", {
  outfile <- tempfile()
  # tmpl_grt_full_path contains the GermlineRootTree marker; see note above
  expect_warning(
    path <- write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
      include_germline_as_tip = TRUE),
    "GermlineRootTree"
  )
  content <- .read_all(path)

  expect_true(grepl(
    '<distribution id="obs.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree.t:run_1">',
    content, fixed = TRUE))
  expect_true(grepl('<taxon id="seq1" spec="Taxon"/>', content, fixed = TRUE))
  expect_true(grepl(
    '<distribution id="germ1.prior" spec="beast.base.evolution.tree.MRCAPrior" tipsonly="true" tree="@Tree.t:run_1">',
    content, fixed = TRUE))
  expect_true(grepl('<taxon id="Germline" spec="Taxon"/>', content, fixed = TRUE))
  expect_true(grepl('lower = "-10000" upper="10000"', content, fixed = TRUE))
})

test_that("${MRCA}: start_date adds a height prior", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, start_date = 500)
  content <- .read_all(path)
  expect_true(grepl('id="height.prior"', content, fixed = TRUE))
  expect_true(grepl('name="mu">500</parameter>', content, fixed = TRUE))
})

test_that("${MRCA}: max_start_date adds a maximum height prior", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, max_start_date = 1000)
  content <- .read_all(path)
  expect_true(grepl('lower="1000" upper="Infinity"', content, fixed = TRUE))
})

# ---- ${ROOTFREQS} -------------------------------------------------------------

test_that("${ROOTFREQS}: empty and no TreeLikelihood rewrite when include_germline_as_root=FALSE", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, include_germline_as_root = FALSE)
  content <- .read_all(path)
  expect_true(grepl('spec="TreeLikelihood"', content, fixed = TRUE))
  expect_false(grepl("rootfreqs.TreeLikelihood", content, fixed = TRUE))
  expect_false(grepl("rootfreqseq", content, fixed = TRUE))
})

test_that("${ROOTFREQS}: rewrites TreeLikelihood -> rootfreqs.TreeLikelihood and embeds germline (no ambiguous bases)", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, include_germline_as_root = TRUE)
  content <- .read_all(path)
  expect_true(grepl('spec="rootfreqs.TreeLikelihood"', content, fixed = TRUE))
  expect_false(grepl('spec="TreeLikelihood"', content, fixed = TRUE)) # fully rewritten, none left
  expect_true(grepl('taxon="Root"', content, fixed = TRUE))
  expect_true(grepl('value="AACCGGTT"/>', content, fixed = TRUE))
})

test_that("${ROOTFREQS}: rootfreqs_class sets the class that replaces TreeLikelihood, rest of the XML is unchanged", {
  write_with <- function(...) {
    path <- write_clone_to_xml(clone1, file = tempfile(), id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
      include_germline_as_tip = FALSE, include_germline_as_root = TRUE, ...)
    .read_all(path)
  }
  tyche_class <- "tyche.evolution.likelihood.RootFreqsTreeLikelihood"
  old <- write_with()                                   # default: today's behavior
  explicit_old <- write_with(rootfreqs_class = "rootfreqs.TreeLikelihood")
  new <- write_with(rootfreqs_class = tyche_class)

  expect_true(grepl('spec="rootfreqs.TreeLikelihood"', old, fixed = TRUE))
  expect_identical(old, explicit_old)
  expect_true(grepl(paste0('spec="', tyche_class, '"'), new, fixed = TRUE))
  expect_false(grepl("rootfreqs.TreeLikelihood", new, fixed = TRUE))
  expect_false(grepl('spec="TreeLikelihood"', new, fixed = TRUE))

  # the <rootfreqseq> element and everything else is the same: only the spec differs
  expect_identical(
    gsub(tyche_class, "rootfreqs.TreeLikelihood", new, fixed = TRUE),
    old)
  expect_true(grepl('<rootfreqseq id="seq_Root" spec="Sequence" taxon="Root"', new, fixed = TRUE))
  expect_true(grepl('value="AACCGGTT"/>', new, fixed = TRUE))
})

test_that("${ROOTFREQS}: rootfreqs_class is ignored when include_germline_as_root=FALSE, and bad names are rejected", {
  path <- write_clone_to_xml(clone1, file = tempfile(), id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, include_germline_as_root = FALSE,
    rootfreqs_class = "tyche.evolution.likelihood.RootFreqsTreeLikelihood")
  content <- .read_all(path)
  expect_true(grepl('spec="TreeLikelihood"', content, fixed = TRUE))
  expect_false(grepl("RootFreqsTreeLikelihood", content, fixed = TRUE))

  for (bad in list("", NA_character_, c("a.B", "c.D"), 1, 'x"/><evil', "a\\1b")) {
    expect_error(
      write_clone_to_xml(clone1, file = tempfile(), id = "run",
        template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
        time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
        include_germline_as_tip = FALSE, include_germline_as_root = TRUE,
        rootfreqs_class = bad),
      "rootfreqs_class must be a single BEAST class name")
  }
})

test_that("${ROOTFREQS}: germline containing 'N' is expanded into per-site probability vectors", {
  clone_n <- .new_toy_clone(clone_id = "2", germline = "AANCCGGT")
  outfile <- tempfile()
  path <- write_clone_to_xml(clone_n, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, include_germline_as_root = TRUE)
  content <- .read_all(path)

  value_str <- .extract(content, 'id="seq_Rootrun_2" spec="Sequence" taxon="Rootrun_2" uncertain="true"\\s*totalcount="4" value="[^"]*"')
  expect_false(is.na(value_str))
  expect_false(grepl("N", value_str, fixed = TRUE)) # 'N' fully replaced

  # exactly one ambiguous site -> exactly one "0.25,0.25,0.25,0.25" group
  n_ambig <- .count_matches(value_str, "0\\.25,0\\.25,0\\.25,0\\.25")
  expect_equal(n_ambig, 1)

  # 8 germline positions -> 8 comma-groups -> 7 semicolons remain after the
  # trailing one is stripped by substr(freqs, 1, nchar(freqs)-1)
  n_semicolons <- .count_matches(value_str, ";")
  expect_equal(n_semicolons, 7)
})

# ---- ${EMP_EQFREQS} -----------------------------------------------------------

test_that("${EMP_EQFREQS}: empirical A/C/G/T frequencies computed from clone@data$sequence", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE)
  content <- .read_all(path)

  # clone1@data$sequence = "AAAAAAAA","CCCCCCCC","GGGGTTTT" -> A=8,C=8,G=4,T=4 of 24
  counts <- c(A = 8, C = 8, G = 4, T = 4)
  expected_freqs <- counts / sum(counts)
  expected_freqstring <- paste(expected_freqs, collapse = " ")

  expect_true(grepl(expected_freqstring, content, fixed = TRUE))
})

# ---- ${ROOT_TYPE_PROBABILITIES} -----------------------------------------------

test_that("${ROOT_TYPE_PROBABILITIES}: character root_trait produces a one-hot vector at the right index", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, trait_list = trait_list_toy, root_trait = "MBC")
  content <- .read_all(path)
  # trait_list_toy = c("GC","MBC") -> MBC is index 2 -> one-hot "0 1"
  expect_true(grepl('name="frequencies">0 1</parameter>', content, fixed = TRUE))
})

test_that("${ROOT_TYPE_PROBABILITIES}: numeric root_trait is treated as a 0-based index", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, trait_list = trait_list_toy, root_trait = 0)
  content <- .read_all(path)
  expect_true(grepl('name="frequencies">1 0</parameter>', content, fixed = TRUE))
})

test_that("${ROOT_TYPE_PROBABILITIES}: dimension reflects length(trait_list), not hardcoded to 2", {
  # trait_data_type_toy/trait passed here are still the 2-state toy trait --
  # trait_list can be fixed to a 3-state list
  outfile <- tempfile()
  trait_list_3state <- c("GC", "MBC", "Plasma")
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, trait_list = trait_list_3state, root_trait = "Plasma")
  content <- .read_all(path)
  # Plasma is index 3 of 3 -> one-hot "0 0 1", dimension must be 3, not 2
  expect_true(grepl('dimension="3" name="frequencies">0 0 1</parameter>', content, fixed = TRUE))
  expect_false(grepl('dimension="2" name="frequencies">0 0 1</parameter>', content, fixed = TRUE))
})


test_that("${ROOT_TYPE_PROBABILITIES}: root_trait not in trait_list raises an error", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
      time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
      include_germline_as_tip = FALSE, trait_list = trait_list_toy, root_trait = "XYZ"),
    "not found in trait_list"
  )
})

test_that("${ROOT_TYPE_PROBABILITIES}: root_trait=NULL strips the placeholder to an empty string", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, root_trait = NULL)
  content <- .read_all(path)
  expect_false(grepl("rootFrequencies", content, fixed = TRUE))
})

test_that("root_trait given but placeholder absent from template raises only a warning", {
  outfile <- tempfile()
  expect_warning(
    path <- write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = FALSE, root_trait = "MBC"),
    "ROOT_TYPE_PROBABILITIES.*not found"
  )
  expect_true(file.exists(path))
})

# ---- ${NODES_TYPE_INIT} -------------------------------------------------------

test_that("${NODES_TYPE_INIT}: defaults to 0 when not supplied via kwargs", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE)
  content <- .read_all(path)
  expect_true(grepl('id="nodeTypesInit" value="0"', content, fixed = TRUE))
})

test_that("${NODES_TYPE_INIT}: uses the value passed via kwargs when supplied", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_full_path, mcmc_length = 5000, log_every = 100,
    time = "collection_date", trait = "celltype", trait_data_type = trait_data_type_toy,
    include_germline_as_tip = FALSE, NODES_TYPE_INIT = "3")
  content <- .read_all(path)
  expect_true(grepl('id="nodeTypesInit" value="3"', content, fixed = TRUE))
})

# ---- <init> / starting tree replacement --------------------------------------

test_that("Providing starting trees: <init> tag is not required at all when no tree is provided", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_broken_no_init_path, mcmc_length = 5000, log_every = 100,
    include_germline_as_tip = FALSE)
  expect_true(file.exists(path))
})

test_that("Providing starting trees: missing <init> tag raises an error only when a starting tree IS provided", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_broken_no_init_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = FALSE, tree = "(seq1:1,seq2:1,seq3:1);"),
    "Could not find <init> tag"
  )
})

test_that("Providing starting trees: a newick character string is used verbatim as the starting tree", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
    include_germline_as_tip = FALSE, tree = "(seq1:1,seq2:1,seq3:1);")
  content <- .read_all(path)
  expect_false(grepl("beast.base.evolution.tree.coalescent.RandomTree", content, fixed = TRUE))
  expect_true(grepl('spec="beast.base.evolution.tree.TreeParser"', content, fixed = TRUE))
  expect_true(grepl('id="NewickTree.t:run_1"', content, fixed = TRUE))
  expect_true(grepl('newick="(seq1:1,seq2:1,seq3:1);"', content, fixed = TRUE))
})

test_that("Providing starting trees: a phylo starting tree drops the Germline tip when include_germline_as_tip=FALSE", {
  toy_tree <- ape::read.tree(text = "(seq1:1,(seq2:1,seq3:1):1,Germline:1);")
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
    include_germline_as_tip = FALSE, tree = toy_tree, tree_states = FALSE)
  content <- .read_all(path)
  newick_attr <- .extract(content, 'newick="[^"]*"')
  expect_false(is.na(newick_attr))
  expect_false(grepl("Germline", newick_attr, fixed = TRUE))
  expect_true(grepl("seq1", newick_attr, fixed = TRUE))
})

test_that("Providing starting trees: a phylo starting tree keeps the Germline tip when include_germline_as_tip=TRUE", {
  toy_tree <- ape::read.tree(text = "(seq1:1,(seq2:1,seq3:1):1,Germline:1);")
  outfile <- tempfile()
  # tmpl_grt_minimal_path contains the GermlineRootTree marker; see note above
    expect_warning(
        path <- write_clone_to_xml(clone1, file = outfile, id = "run",
        template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
        include_germline_as_tip = TRUE, tree = toy_tree, tree_states = FALSE),
    "GermlineRootTree"
  )
  content <- .read_all(path)
  newick_attr <- .extract(content, 'newick="[^"]*"')
  expect_false(is.na(newick_attr))
  expect_true(grepl("Germline", newick_attr, fixed = TRUE))
})

test_that("Providing starting trees: a phylo starting tree correctly creates a Newick tree", {
  # this is an actually rooted tree, not just a trifurcation, to test that the correct tree is preserved in the output
  toy_tree <- ape::read.tree(text = "(((seq2:1,seq3:1):1,seq1:1):1,Germline:1);")
  outfile <- tempfile()
  # tmpl_grt_minimal_path contains the GermlineRootTree marker; see note above
    expect_warning(
        path <- write_clone_to_xml(clone1, file = outfile, id = "run",
        template = tmpl_grt_minimal_path, mcmc_length = 5000, log_every = 100,
        include_germline_as_tip = TRUE, tree = toy_tree, tree_states = FALSE, start_edge_length = 1),
    "GermlineRootTree"
  )
  content <- .read_all(path)
  newick_attr <- .extract(content, 'newick="[^"]*"')
  expect_false(is.na(newick_attr))
  output_tree <- ape::read.tree(text=newick_attr)
  expect_true(ape::all.equal.phylo(toy_tree, output_tree))
})

test_that("Providing starting trees: tree_states=TRUE back-fills nodeTypes with zeros when tree has no $state", {
  toy_tree <- ape::read.tree(text = "(seq1:1,(seq2:1,seq3:1):1,Germline:1);")
  outfile <- tempfile()
  expect_warning(
    path <- write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_minimal_tc_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = TRUE, tree = toy_tree, tree_states = TRUE),
    "States not found"
  )
  content <- .read_all(path)
  # nrow(clone1@data) == 3 -> tips = 4 -> nodes = 7 -> 7 zeros
  expect_true(grepl("id='nodeTypes'", content, fixed = TRUE))
  expect_true(grepl('value="0 0 0 0 0 0 0"', content, fixed = TRUE))
})

test_that("Providing starting trees: tree_states=TRUE errors if no nodeTypes parameter exists in the template", {
  toy_tree <- ape::read.tree(text = "(seq1:1,(seq2:1,seq3:1):1,Germline:1);")
  outfile <- tempfile()
  expect_error(
    suppressWarnings(write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_grt_minimal_no_tc_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = TRUE, tree = toy_tree, tree_states = TRUE)),
    "Could not find nodeTypes"
  )
})

# ---- generic catch-all ${VAR} + deferred `replacements` -----------------------

test_that("free-form placeholders filled from kwargs", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_catchall_path, mcmc_length = 5000, log_every = 100,
    include_germline_as_tip = FALSE, FOO = "hello", BAR = "world")
  content <- .read_all(path)
  expect_true(grepl('value="hello"', content, fixed = TRUE))
  expect_true(grepl('value="world"', content, fixed = TRUE))
})

test_that("an unfilled placeholder raises an error naming that variable", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_catchall_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = FALSE, FOO = "hello"),
    "BAR"
  )
})

test_that("`replacements`: deferred variable is substituted after the main loop, verbatim (not re-scanned)", {
  outfile <- tempfile()
  path <- write_clone_to_xml(clone1, file = outfile, id = "run",
    template = tmpl_catchall_path, mcmc_length = 5000, log_every = 100,
    include_germline_as_tip = FALSE, replacements = "BAR",
    BAR = "ref-${FOO}-end", FOO = "hello")
  content <- .read_all(path)
  # FOO substituted normally in the main loop
  expect_true(grepl('value="hello"', content, fixed = TRUE))
  # BAR substituted afterward, verbatim -- the literal "${FOO}" text
  # embedded in BAR's own value is NOT re-processed a second time
  expect_true(grepl('value="ref-${FOO}-end"', content, fixed = TRUE))
})

# ---- write_clones_to_xmls(): multi-clone wrapper behavior ---------------------

test_that("write_clones_to_xmls(): trait_list is the union across all clones, shared codeMap in every file", {
  clone_x <- .new_toy_clone(clone_id = "10", germline = "AACCGGTT") # celltype: GC,GC,MBC
  clone_y <- .new_toy_clone(clone_id = "20", germline = "AACCGGTT")
  clone_y@data$celltype <- c("GC", "Plasma", "Plasma")

  outfile <- tempfile()
  paths <- write_clones_to_xmls(list(clone_x, clone_y), id = "wrap",
    outfile = outfile, template = tmpl_wrap_path, trait = "celltype",
    include_germline_as_tip = FALSE, mcmc_length = 5000, log_every = 100)

  expect_equal(paths, c(paste0(outfile, "_10.xml"), paste0(outfile, "_20.xml")))
  expect_true(all(file.exists(paths)))

  expected_trait_list <- sort(unique(c("GC", "MBC", "Plasma")))
  expected_data_type <- .build_trait_data_type(expected_trait_list)

  content_x <- .read_all(paths[1])
  content_y <- .read_all(paths[2])

  # same codeMap embedded identically in BOTH files, computed once from the
  # union of trait values across clones -- not recomputed per-clone
  expect_true(grepl(expected_data_type, content_x, fixed = TRUE))
  expect_true(grepl(expected_data_type, content_y, fixed = TRUE))

  # each file still reflects its own clone's actual trait values
  expect_true(grepl("seq3=MBC", content_x, fixed = TRUE))
  expect_true(grepl("seq3=Plasma", content_y, fixed = TRUE))
})

test_that("write_clones_to_xmls(): each clone gets its own starting tree from the `trees` list, by position", {
  clone_x <- .new_toy_clone(clone_id = "10", germline = "AACCGGTT")
  clone_y <- .new_toy_clone(clone_id = "20", germline = "AACCGGTT")
  clone_y@data$celltype <- c("GC", "Plasma", "Plasma")

  tree_x <- "(seq1:1,seq2:1,seq3:1);"
  tree_y <- "(seq1:2,seq2:2,seq3:2);"

  outfile <- tempfile()
  paths <- write_clones_to_xmls(list(clone_x, clone_y), id = "wrap",
    outfile = outfile, template = tmpl_wrap_path, trait = "celltype",
    include_germline_as_tip = FALSE, mcmc_length = 5000, log_every = 100,
    trees = list(tree_x, tree_y))

  content_x <- .read_all(paths[1])
  content_y <- .read_all(paths[2])

  expect_true(grepl(paste0('newick="', tree_x, '"'), content_x, fixed = TRUE))
  expect_true(grepl(paste0('newick="', tree_y, '"'), content_y, fixed = TRUE))
  # not swapped
  expect_false(grepl(paste0('newick="', tree_y, '"'), content_x, fixed = TRUE))
  expect_false(grepl(paste0('newick="', tree_x, '"'), content_y, fixed = TRUE))
})

# ============================================================================================
# Part 2: describeTemplate() and the dowser-template-info convention
# ============================================================================================
# The helpers .write_template(), .new_toy_clone() and clone1 are the ones at the top of this file.

# ---- Toy templates ----------------------------------------------------------
#
# DOC: every structural element the built-in checklist recognizes is
#      present, plus a dowser-template-info block documenting two free-form
#      vars (one with a wrapped, multi-line description) and one @struct
#      note.
# UNDOC: the structural opposite of DOC (every checklist element absent
#        except <init>, which every valid template needs), no
#        dowser-template-info block at all, one undocumented free-form var.
# MALFORMED: has the opening marker but never closes with '-->'.
# HOOK: minimal catch-all template with FOO documented and BAR left
#       undocumented, for testing write_clone_to_xml()'s error-hook.
# NO_DATA: missing ${DATA} entirely, for testing the var="DATA" lookup's
#          absent branch -- not a usable template (every real template
#          needs ${DATA}), but that's exactly the case this branch reports.

tmpl_doc_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<parameter id="nodeTypes" dimension="${NODES}">0</parameter>',
  '<state>',
  '<tree id="Tree.t:${CLONE}" spec="tyche.evolution.tree.GermlineRootTree" name="stateNode">',
  '${DATE}',
  '</tree>',
  '<parameter id=\'nodeTypes\' name="stateNode">0</parameter>',
  '</state>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '<distribution id="prior" spec="CompoundDistribution">',
  '${MRCA}',
  '<distribution id="typedTreeLikelihood.type" spec="tyche.evolution.likelihood.AncestralTypeLikelihood" tag="${TRAIT_NAME}" tree="@Tree.t:${CLONE}">',
  '${TRAIT}',
  '${ROOT_TYPE_PROBABILITIES}',
  '</distribution>',
  '</distribution>',
  '<distribution id="treeLikelihood.${CLONE}" spec="TreeLikelihood" data="@${CLONE}" tree="@Tree.t:${CLONE}">',
  '${ROOTFREQS}',
  '<parameter id="freqParameter.s:${CLONE}" dimension="4">${EMP_EQFREQS}</parameter>',
  '</distribution>',
  '<parameter id="nodeTypesInit" value="${NODES_TYPE_INIT}"/>',
  '<parameter id="fooRateParam" value="${FOO_RATE}"/>',
  '<parameter id="barRateParam" value="${BAR_RATE}"/>',
  '${OPERATORS}',
  '</run>',
  '<!-- dowser-template-info',
  '@var FOO_RATE | required | Made-up rate parameter for testing, part one',
  '  of a description that wraps onto a second line.',
  '@var BAR_RATE | optional, default 1.0 | Made-up rate parameter for testing.',
  '@struct GermlineRootTree | Test author note on GRT.',
  '-->',
  '</beast>'
))

tmpl_undoc_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<state>',
  '<tree id="Tree.t:${CLONE}" spec="beast.base.evolution.tree.Tree" name="stateNode">',
  '</tree>',
  '</state>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '<parameter id="bazRateParam" value="${BAZ_RATE}"/>',
  '</run>',
  '</beast>'
))

tmpl_malformed_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '</run>',
  '<!-- dowser-template-info',
  '@var FOO | required | this block never closes',
  '</beast>'
))

tmpl_hook_path <- .write_template(c(
  '<beast>',
  '${DATA}',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<parameter id="fooParam.${CLONE}" value="${FOO}"/>',
  '<parameter id="barParam.${CLONE}" value="${BAR}"/>',
  '<init id="RandomTree.t:${CLONE}" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:${CLONE}" taxa="@${CLONE}">',
  '</init>',
  '<logger id="tracelog" fileName="${CLONE}.log" logEvery="1000">',
  '</logger>',
  '<logger id="treeWithTraitLogger" fileName="${CLONE}.trees" logEvery="1000" mode="tree">',
  '</logger>',
  '</run>',
  '<!-- dowser-template-info',
  '@var FOO | required | The foo rate, for testing the error-hook lookup.',
  '-->',
  '</beast>'
))

tmpl_no_data_path <- .write_template(c(
  '<beast>',
  '<run id="mcmc" spec="MCMC" chainLength="30000000">',
  '<parameter id="fooParam.${CLONE}" value="${FOO}"/>',
  '</run>',
  '</beast>'
))

# ---- .parseTemplateInfo() ---------------------------------------------------

test_that(".parseTemplateInfo(): parses @var entries, joining a wrapped continuation line", {
  info <- .parseTemplateInfo(readLines(tmpl_doc_path))

  expect_true(info$found)
  expect_equal(info$vars$FOO_RATE$tier, "required")
  expect_equal(
    info$vars$FOO_RATE$description,
    "Made-up rate parameter for testing, part one of a description that wraps onto a second line."
  )
  expect_equal(info$vars$BAR_RATE$tier, "optional, default 1.0")
  expect_equal(info$vars$BAR_RATE$description, "Made-up rate parameter for testing.")
})

test_that(".parseTemplateInfo(): parses @struct notes", {
  info <- .parseTemplateInfo(readLines(tmpl_doc_path))
  expect_equal(info$structs$GermlineRootTree, "Test author note on GRT.")
})

test_that(".parseTemplateInfo(): found=FALSE and empty lists when no block exists", {
  info <- .parseTemplateInfo(readLines(tmpl_undoc_path))
  expect_false(info$found)
  expect_equal(length(info$vars), 0)
  expect_equal(length(info$structs), 0)
})

test_that(".parseTemplateInfo(): warns and ignores a block with no closing '-->'", {
  expect_warning(
    info <- .parseTemplateInfo(readLines(tmpl_malformed_path)),
    "no closing"
  )
  expect_false(info$found)
})

# ---- .templateVarLookup() ---------------------------------------------------

test_that(".templateVarLookup(): returns the documented description when one exists", {
  hint <- .templateVarLookup(readLines(tmpl_doc_path), "FOO_RATE")
  expect_true(grepl("Made-up rate parameter", hint, fixed = TRUE))
})

test_that(".templateVarLookup(): returns NULL for an undocumented variable", {
  hint <- .templateVarLookup(readLines(tmpl_undoc_path), "BAZ_RATE")
  expect_null(hint)
})

test_that(".templateVarLookup(): returns NULL when there's no dowser-template-info block at all", {
  hint <- .templateVarLookup(readLines(tmpl_undoc_path), "anything")
  expect_null(hint)
})

# ---- describeTemplate(): structural checklist -------------------------------

test_that("describeTemplate(): every built-in structural element is [present] when its tag/placeholder exists", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE)
  expect_true(all(result$structural$present))
})

test_that("describeTemplate(): every built-in structural element is [absent] except <init>, which every valid template needs", {
  result <- describeTemplate(tmpl_undoc_path, print = FALSE)
  init_row <- result$structural[result$structural$element == "<init> tag", ]
  other_rows <- result$structural[result$structural$element != "<init> tag", ]

  expect_true(init_row$present)
  expect_true(all(!other_rows$present))
})

test_that("describeTemplate(): an author's @struct note is appended to the matching built-in row", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE)
  grt_row <- result$structural[grepl("GermlineRootTree", result$structural$element, fixed = TRUE), ]
  expect_true(grepl("Author note: Test author note on GRT.", grt_row$effect, fixed = TRUE))
})

# ---- describeTemplate(): variables ------------------------------------------

test_that("describeTemplate(): documented free-form vars are reported with their tier and description", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE)
  foo_row <- result$vars[result$vars$name == "FOO_RATE", ]
  bar_row <- result$vars[result$vars$name == "BAR_RATE", ]

  expect_true(foo_row$documented)
  expect_equal(foo_row$tier, "required")
  expect_true(bar_row$documented)
  expect_equal(bar_row$tier, "optional, default 1.0")
})

test_that("describeTemplate(): an undocumented free-form var is flagged, not silently dropped", {
  result <- describeTemplate(tmpl_undoc_path, print = FALSE)
  baz_row <- result$vars[result$vars$name == "BAZ_RATE", ]

  expect_equal(nrow(baz_row), 1)
  expect_false(baz_row$documented)
  expect_equal(baz_row$tier, "undocumented")
  expect_equal(baz_row$description, "")
})

test_that("describeTemplate(): Dowser-special placeholders (DATA, CLONE, DATE, etc.) are excluded from the vars table", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE)
  expect_false(any(result$vars$name %in% c(
    "DATA", "CLONE", "DATE", "TRAIT", "TRAIT_NAME", "NODES", "MRCA",
    "OPERATORS", "ROOTFREQS", "EMP_EQFREQS", "ROOT_TYPE_PROBABILITIES",
    "NODES_TYPE_INIT"
  )))
})

# ---- describeTemplate(): printing behavior and errors -----------------------

test_that("describeTemplate(): print=FALSE produces no console output", {
  output <- capture.output(describeTemplate(tmpl_doc_path, print = FALSE))
  expect_equal(length(output), 0)
})

test_that("describeTemplate(): print=TRUE (default) prints present/absent markers and var names", {
  output <- paste(capture.output(describeTemplate(tmpl_doc_path)), collapse = "\n")
  expect_true(grepl("[present]", output, fixed = TRUE))
  expect_true(grepl("FOO_RATE", output, fixed = TRUE))
  expect_true(grepl("All 2 variables are documented", output, fixed = TRUE))
})

test_that("describeTemplate(): print=TRUE flags undocumented variables in the printed summary", {
  output <- paste(capture.output(describeTemplate(tmpl_undoc_path)), collapse = "\n")
  expect_true(grepl("[absent]", output, fixed = TRUE))
  expect_true(grepl("UNDOCUMENTED", output, fixed = TRUE))
  expect_true(grepl("1 of 1 variables are undocumented", output, fixed = TRUE))
})

test_that("describeTemplate(): errors clearly if the file doesn't exist", {
  expect_error(
    describeTemplate(tempfile(fileext = ".xml"), print = FALSE),
    "does not exist"
  )
})

# ---- describeTemplate(..., var = ): single-variable lookup ------------------

test_that("describeTemplate(var=): a documented free-form variable reports its tier and description", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE, var = "FOO_RATE")
  expect_equal(result$kind, "free")
  expect_true(result$present)
  expect_true(result$documented)
  expect_equal(result$tier, "required")
  expect_true(grepl("Made-up rate parameter", result$description, fixed = TRUE))
})

test_that("describeTemplate(var=): an undocumented free-form variable is flagged as such, not silently missing", {
  result <- describeTemplate(tmpl_undoc_path, print = FALSE, var = "BAZ_RATE")
  expect_equal(result$kind, "free")
  expect_true(result$present)
  expect_false(result$documented)
  expect_equal(result$tier, "undocumented")
})

test_that("describeTemplate(var=): a present Dowser-special placeholder reports kind='structural' and present=TRUE", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE, var = "EMP_EQFREQS")
  expect_equal(result$kind, "structural")
  expect_true(result$present)
  expect_true(nchar(result$effect) > 0)
})

test_that("describeTemplate(var=): an absent Dowser-special placeholder reports present=FALSE", {
  result <- describeTemplate(tmpl_undoc_path, print = FALSE, var = "MRCA")
  expect_equal(result$kind, "structural")
  expect_false(result$present)
})

test_that("describeTemplate(var=): TRAIT and TRAIT_NAME map to the same combined structural check", {
  trait_result <- describeTemplate(tmpl_doc_path, print = FALSE, var = "TRAIT")
  trait_name_result <- describeTemplate(tmpl_doc_path, print = FALSE, var = "TRAIT_NAME")
  expect_equal(trait_result$present, trait_name_result$present)
  expect_equal(trait_result$effect, trait_name_result$effect)
})

test_that("describeTemplate(var=): DATA reports present=TRUE when ${DATA} is in the template", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE, var = "DATA")
  expect_equal(result$kind, "structural")
  expect_true(result$present)
})

test_that("describeTemplate(var=): DATA reports present=FALSE, with a clear message, when ${DATA} is missing", {
  result <- describeTemplate(tmpl_no_data_path, print = FALSE, var = "DATA")
  expect_equal(result$kind, "structural")
  expect_false(result$present)
  expect_true(grepl("cannot be used by Dowser", result$effect, fixed = TRUE))
})

test_that("describeTemplate(var=): CLONE reports present=TRUE when ${CLONE} is in the template", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE, var = "CLONE")
  expect_equal(result$kind, "structural")
  expect_true(result$present)
})

test_that("describeTemplate(var=): a variable with no ${VAR} anywhere in the file is reported as not_found", {
  result <- describeTemplate(tmpl_doc_path, print = FALSE, var = "TOTALLY_MADE_UP_VAR")
  expect_equal(result$kind, "not_found")
  expect_false(result$present)
})

test_that("describeTemplate(var=): print=FALSE produces no console output", {
  output <- capture.output(describeTemplate(tmpl_doc_path, print = FALSE, var = "FOO_RATE"))
  expect_equal(length(output), 0)
})

test_that("describeTemplate(var=): print=TRUE shows the variable name and template file in the header line", {
  output <- paste(capture.output(describeTemplate(tmpl_doc_path, var = "FOO_RATE")), collapse = "\n")
  expect_true(grepl("FOO_RATE", output, fixed = TRUE))
  expect_true(grepl(basename(tmpl_doc_path), output, fixed = TRUE))
})

test_that("describeTemplate(var=): print=TRUE for a structural variable notes it isn't a free-form variable", {
  output <- paste(capture.output(describeTemplate(tmpl_doc_path, var = "EMP_EQFREQS")), collapse = "\n")
  expect_true(grepl("[present]", output, fixed = TRUE))
  expect_true(grepl("Dowser-recognized, not a free-form variable", output, fixed = TRUE))
})

test_that("describeTemplate(var=): print=TRUE for a not-found variable says so plainly", {
  output <- paste(capture.output(describeTemplate(tmpl_doc_path, var = "TOTALLY_MADE_UP_VAR")), collapse = "\n")
  expect_true(grepl("Not found", output, fixed = TRUE))
})

# ---- error-hook: write_clone_to_xml() uses the doc block in its error -------

test_that("write_clone_to_xml(): missing a documented variable includes its description in the error", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_hook_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = FALSE, BAR = "world"),
    "FOO.*not provided.*The foo rate, for testing the error-hook lookup"
  )
})

test_that("write_clone_to_xml(): missing an undocumented variable points to describeTemplate() instead", {
  outfile <- tempfile()
  expect_error(
    write_clone_to_xml(clone1, file = outfile, id = "run",
      template = tmpl_hook_path, mcmc_length = 5000, log_every = 100,
      include_germline_as_tip = FALSE, FOO = "hello"),
    "BAR.*not provided.*describeTemplate"
  )
})
# nolint end
