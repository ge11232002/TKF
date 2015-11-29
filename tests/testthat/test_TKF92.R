test_that("test_TKF92", {
    library(seqinr)
    fasta <- read.fasta(file.path(system.file("extdata", package = "TKF"), 
        "pair1.fasta"), seqtype = "AA", set.attributes = FALSE)
    seq1 <- fasta[[1]]
    seq2 <- fasta[[2]]
    ans <- TKF92Pair(seq1, seq2, mu = 0.0006137344, r = 0.7016089061, 
        substModel = GONNET, substModelBF = GONNETBF)
    expect_equal(ans["PAM"], 116.6130887028, tolerance=1e-3)
})

