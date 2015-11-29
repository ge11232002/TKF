test_that("test_TKF91", {
    library(seqinr)
    fasta <- read.fasta(file.path(system.file("extdata", package = "TKF"), 
        "pair1.fasta"), seqtype = "AA", set.attributes = FALSE)
    seq1 <- fasta[[1]]
    seq2 <- fasta[[2]]
    ans <- TKF91Pair(seq1, seq2, mu = 0.0005920655, substModel = GONNET, 
        substModelBF = GONNETBF)
    expectedPam <- 116.3416794647
    expect_equal(ans["PAM"], expectedPam, tolerance=1e-3, 
                 check.attributes=FALSE)
})

