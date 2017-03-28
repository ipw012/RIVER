library(RIVER)
context("Two main functions for RIVER")

test_that("Two main functions work as expected", {
  dataInput <-
    getData(filename=system.file("extdata", "simulation_RIVER.gz",
                                 package = "RIVER"), ZscoreThrd=1.5)

  theta_init <- matrix(c(.99, .01, .3, .7), nrow=2)
  costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)

  ## evaRIVER works as expected
  evaROC <- evaRIVER(dataInput)
  expect_equal(length(evaROC$RIVER_sens),length(evaROC$RIVER_spec))
  expect_equal(length(evaROC$GAM_sens),length(evaROC$GAM_spec))
  expect_lte(evaROC$pvalue,1)
  expect_equal(class(evaROC$RIVER_auc),"numeric")

  ## appRIVER works as expected
  postprobs <- appRIVER(dataInput, pseudoc=50, theta_init, costs)
  expect_is(postprobs$indiv_name,"character")
  expect_is(postprobs$gene_name,"character")
  expect_equal(length(postprobs$RIVER_posterior),
               length(postprobs$GAM_posterior))
  expect_identical(sampleNames(dataInput),
                   paste(postprobs$indiv_name,":",
                         postprobs$gene_name,sep=""))
  expect_is(postprobs$fitRIVER, "list")
})

context("Two optional plot functions for RIVER")

test_that("Two plot functions work as expected", {
  dataInput <-
    getData(filename=system.file("extdata", "simulation_RIVER.gz",
                                 package = "RIVER"), ZscoreThrd=1.5)

  theta_init <- matrix(c(.99, .01, .3, .7), nrow=2)
  costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)

  rocstat <- evaRIVER(dataInput, pseudoc=50, theta_init, costs)

  outRIVER <- appRIVER(dataInput, pseudoc=50, theta_init, costs)

  ## plotPosteriors works as expected
  expect_silent(
    plotPosteriors(outRIVER, as.numeric(unlist(dataInput$Outlier))-1))
})
