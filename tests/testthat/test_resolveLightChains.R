# Load test dataset
#library(testthat)
#library(dowser)

data <- data.frame(sequence_id = c(unlist(lapply(1:30, function(x) paste0("seq", x)))),
                   subject_id = "S1",
                   v_call = c("IGHV1-1*01","IGKV1-1*01","IGHV1-1*01","IGKV1-1*01",
                              "IGHV1-1*01","IGKV2-1*01","IGHV1-1*01","IGKV1-1*01",
                              "IGKV2-1*01","IGHV1-1*01","IGHV1-1*01","IGHV1-1*01",
                              "IGKV1-1*01,IGKV2-1*01","IGHV1-1*01","IGKV3-1*01",
                              "IGHV2-1*01","IGKV1-1*01","IGKV1-1*01","IGKV1-1*01",
                              "IGHV1-1*01","IGHV1-1*01","IGHV1-1*01","IGKV4-1*01",
                              "IGHV1-1*01","IGKV4-1*01","IGHV1-1*01","IGHV1-1*01",
                              "IGKV1-1*01","IGHV1-1*01","IGKV1-1*01"),
                   j_call = c("IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01","IGHJ2*01",
                              "IGKJ1*01","IGHJ2*01","IGKJ1*01","IGKJ1*01","IGHJ2*01",
                              "IGHJ2*01","IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01",
                              "IGHJ1*01","IGKJ1*01","IGKJ1*01","IGKJ1*01","IGHJ2*01",
                              "IGHJ2*01","IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01",
                              "IGHJ2*01","IGHJ2*01","IGKJ1*01","IGHJ2*01","IGKJ1*01"),
                   junction = c("TGTAAAAAATGG","TGTCCCCCCTGG","TGTAAAAAATGG","TGTCCCCCCTGG",
                                "TGTAAAAAATGT","TGTCCCCCCTGG","TGTAAAAAATGG","TGTCCCCCCTGG",
                                "TGTCCCCCCTGG","TGTAAAAAATGG","TGTAAAAAATGG","TGTAAAAAATGG",
                                "TGTCCCCCCTGG","TGTAAAAAATCG","TGTCCCCCCTGG","TGTAAAACCTGG",
                                "TGTCCCCCCTGG","TGTCCCCCCTGG","TGTCCCCCCTGG","TGTAAAAAATGT",
                                "TGTAAAAAATCG","TGTAAAAAATGG","TGTCCCCCCTGG","TGTAAAAAATGG",
                                "TGTCCCCCCTGG","TGTAGAAAATGG","TGTAAAAAATGG","TGTCCCCCCTGG",
                                "TGTAAAAAATGG","TGTCCCCCCTGG"),
                   locus = c("IGH","IGK","IGH","IGK","IGH","IGK","IGH","IGK","IGK",
                             "IGH","IGH","IGH","IGK","IGH","IGK","IGH","IGK","IGK",
                             "IGK","IGH","IGH","IGH","IGK","IGH","IGK","IGH","IGH",
                             "IGK","IGH","IGK"),
                   cell_id = c(1,1,2,2,3,3,4,4,4,6,NA,8,8,5,5,7,7,9,NA,10,11,12,
                               12,13,13,14,15,15,16,16),
                   junction_length = 12,
                   clone_id = c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,NA,NA,1,2,4,4,
                                4,4,1,1,1,1,1),
                   expected_clone_subgroup = c(1,1,1,1,2,2,1,1,NA,1,1,1,1,1,1,1,
                                               1,NA,NA,2,1,1,1,1,1,1,1,1,1,1))


expect_warning(data <- resolveLightChains(data, seq = "junction"))

expect_equal("seq9" %in% data$sequence_id, FALSE)
expect_equal("seq18" %in% data$sequence_id, FALSE)
expect_equal("seq19" %in% data$sequence_id, FALSE)
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq5")], data$clone_subgroup[which(data$sequence_id == "seq20")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq15")], data$clone_subgroup[which(data$sequence_id == "seq21")])
expect_equal(data$clone_subgroup[which(data$sequence_id == "seq26")], 1)
