test_that("optim.phylo.wls", {
  distFn <- file.path(system.file("extdata", package="TKF"),
                      "r2.distanceMatrix")
  Dist <- read.table(distFn, header=TRUE, sep="\t")
  rownames(Dist) <- colnames(Dist)
  Dist <- as.matrix(Dist)
  varFn <- file.path(system.file("extdata", package="TKF"),
                     "r2.varianceMatrix")
  Var <- read.table(varFn, header=TRUE, sep="\t")
  rownames(Var) <- colnames(Var)
  Var <- as.matrix(Var)
  tree <- optim.phylo.wls(Dist, Var)
  expect_is(tree, "phylo")
  })