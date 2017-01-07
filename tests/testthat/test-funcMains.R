library(RIVERpkg)
context("Two main functions for RIVER")

test_that("Two main functions work as expected", {
  ## simulated data
  G <- simulated_features
  E <- simulated_outliers

  ## check input data
  expect_identical(G[,"SubjectID"],E[,"SubjectID"])
  expect_identical(G[,"GeneName"],E[,"GeneName"])

  ## experimental setup
  theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
  costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)

  ## check experimental setup
  expect_equal(dim(G)[1],dim(E)[1])
  expect_equal(sum(unique(E[,"Outlier"])),1)
  expect_equal(colSums(theta_init),c(1,1))

  rocstat <- evaRIVER(G, E, pseudoc=50, theta_init, costs)

  ## evaRIVER works as expected
  expect_is(rocstat, "list")
  expect_equal(length(rocstat),3)
  expect_equal(length(rocstat$rocRIVER$specificities),
               length(rocstat$rocRIVER$sensitivities))
  expect_equal(length(rocstat$rocGAM$specificities),
               length(rocstat$rocGAM$sensitivities))
  expect_lte(rocstat$pValue,1)

  outRIVER <- appRIVER(G, E, pseudoc=50, theta_init, costs)

  ## appRIVER works as expected
  expect_is(outRIVER, "list")
  expect_equal(length(outRIVER),5)
  expect_equal(length(outRIVER$pFRgivenG),dim(G)[1])
  expect_equal(length(outRIVER$pFRgivenGE),dim(G)[1])
})

context("Two optional plot functions for RIVER")

test_that("Two plot functions work as expected", {
  ## simulated data
  G <- simulated_features
  E <- simulated_outliers

  ## experimental setup
  theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
  costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)

  rocstat <- evaRIVER(G, E, pseudoc=50, theta_init, costs)
  ## plotAUC works as expected
  expect_silent(plotAUC(rocstat))

  outRIVER <- appRIVER(G, E, pseudoc=50, theta_init, costs)

  ## plotPosteriors works as expected
  expect_silent(plotPosteriors(outRIVER, E))
})
