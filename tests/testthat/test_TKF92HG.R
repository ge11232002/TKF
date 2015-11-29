test_that("test_TKF92HG", {
    library(seqinr)
    fasta <- read.fasta(file.path(system.file("extdata", package = "TKF"), 
        "pair1.fasta"), seqtype = "AA", set.attributes = FALSE)
    seq1 <- fasta[[1]]
    seq2 <- fasta[[2]]
    ans <- TKF92HGPair(seq1, seq2, mu = 0.0005920655, r = 0.8, 
        Ps = 1, Kf = 1.2, substModel = GONNET, substModelBF = GONNETBF)
    expectedPam <- 119.3832517
    expect_equal(ans["PAM"], expectedPam, tolerance=1e-3, 
                 check.attributes=FALSE)
})

