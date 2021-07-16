#library(testthat)
#library(dowser)

test_that("makeAirrClone",{
    # Example AIRR data.frame
    db <- data.frame(sequence_id=LETTERS[1:4],
                     sequence_alignment=c("CCCCTGGG", "CCCCTGGN", "NAACTGGN", "NNNCTGNN"),
                     v_call="Homsap IGKV1-39*01 F",
                     j_call="Homsap IGKJ5*01 F",
                     junction_length=2,
                     germline_alignment="CCCCAGGG",
                     clone_id=1,
                     isotype=c("IgG", "IgG", "IgM", "IgA"),
                     stringsAsFactors=FALSE)

    exp <- data.frame("sequence_id"=c("C", "A"),
                      "sequence"=c("NAACTGGNN", "CCCCTGGGN"),
                      "lsequence"=c("", ""),
                      "hlsequence"=c("NAACTGGNN", "CCCCTGGGN"),
                      "collapse_count"=c(1, 2),
                      stringsAsFactors=FALSE)

    exp_trait <- data.frame("sequence_id"=c("D", "A", "C"),
                      "sequence"=c("NNNCTGNNN","CCCCTGGGN", "NAACTGGNN"),
                      "isotype"=c("IgA","IgG","IgM"),
                      "lsequence"=c("", "", ""),
                      "hlsequence"=c("NNNCTGNNN","CCCCTGGGN", "NAACTGGNN"),
                      "collapse_count"=c(1, 2, 1),
                      stringsAsFactors=FALSE)
    
    # Without splitting by trait value
    clones <- formatClones(db,germ="germline_alignment",randomize=FALSE,
      use_regions=FALSE)
    clone <- clones$data[[1]]
    
    expect_equal(clone@clone, "1")
    expect_equal(clone@germline, "CCCCAGGGN")
    expect_equal(clone@v_gene, "IGKV1-39")
    expect_equal(clone@j_gene, "IGKJ5")
    expect_equal(clone@junc_len, 2)
    expect_equal(clone@locus,rep("N",length=9))
    expect_equal(clone@phylo_seq,"sequence")
    expect_equal(clone@lgermline,"")
    expect_equal(clone@hlgermline,"CCCCAGGGN")
    expect_equal(clone@germline,"CCCCAGGGN")
    expect_equal(clone@data, exp, tolerance=0.001)

    # With splitting by trait value
    clones <- formatClones(db,germ="germline_alignment",randomize=FALSE,
      trait="isotype",add_count=TRUE,use_regions=FALSE)
    clone <- clones$data[[1]]
    
    expect_equal(clone@clone, "1")
    expect_equal(clone@germline, "CCCCAGGGN")
    expect_equal(clone@v_gene, "IGKV1-39")
    expect_equal(clone@j_gene, "IGKJ5")
    expect_equal(clone@junc_len, 2)
    expect_equal(clone@locus,rep("N",length=9))
    expect_equal(clone@phylo_seq,"sequence")
    expect_equal(clone@lgermline,"")
    expect_equal(clone@hlgermline,"CCCCAGGGN")
    expect_equal(clone@germline,"CCCCAGGGN")
    expect_equal(clone@data, exp_trait, tolerance=0.001)
})
