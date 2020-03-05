library(RIVER)
context("Basic functions for RIVER")

test_that("Basic functions work as expected", {
  ## getData
  dataInput <- getData(filename=system.file("extdata", "simulation_RIVER.gz",
                                           package = "RIVER"), ZscoreThrd=1.5)
  FeatAll <- t(exprs(dataInput)) # all genomic features
  OutAll <- as.numeric(unlist(dataInput$Outlier))-1 # all outlier status
  expect_equal(inherits(FeatAll,"matrix"),TRUE)
  expect_equal(inherits(OutAll,"numeric"),TRUE)
  expect_equal(length(unique(OutAll)),2) # either outlier or nonoutlier

  # training G data
  FeatTrng <- t(exprs(dataInput[,is.na(dataInput$N2pair)]))
  # training E data
  OutTrng <-
    as.numeric(unlist(dataInput$Outlier[is.na(dataInput$N2pair)]))-1
  expect_equal(inherits(FeatTrng,"matrix"),TRUE)
  expect_equal(inherits(OutTrng,"numeric"),TRUE)

  # Test G data
  FeatTest <-
    t(cbind(exprs(dataInput[,!is.na(dataInput$N2pair)])
            [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
            exprs(dataInput[,!is.na(dataInput$N2pair)])
            [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))
  # Test E data (indiv1 and then indiv2 from N2 pairs)
  OutTest1 <-
    as.numeric(unlist(
      c(dataInput$Outlier[!is.na(dataInput$N2pair)]
        [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
        dataInput$Outlier[!is.na(dataInput$N2pair)]
        [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1
  # Test E data (indiv2 and then indiv1 from N2 pairs)
  OutTest2 <-
    as.numeric(unlist(
      c(dataInput$Outlier[!is.na(dataInput$N2pair)]
        [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
        dataInput$Outlier[!is.na(dataInput$N2pair)]
        [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1
  expect_equal(inherits(FeatTest,"matrix"),TRUE)
  expect_identical(
    rownames(FeatTest),
    c(sampleNames(dataInput[,!is.na(dataInput$N2pair)])
      [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
      sampleNames(dataInput[,!is.na(dataInput$N2pair)])
      [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))
  expect_equal(length(unique(OutTest1)),2)
  expect_equal(length(unique(OutTest2)),2)

  # N2 pairs should have same genomic features
  expect_equivalent(
    sum(diag(cor(exprs(dataInput[,!is.na(dataInput$N2pair)])
                 [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                 exprs(dataInput[,!is.na(dataInput$N2pair)])
                 [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))),
    sum(!is.na(dataInput$N2pair))/2)

  theta_init <- matrix(c(.99, .01, .3, .7), nrow=2)
  costs <- c(100, 10, 1, .1, .01, 1e-3, 1e-4)

  ## check experimental setup
  expect_equal(sum(unique(OutAll)),1)
  expect_equal(colSums(theta_init),c(1,1))

  logisticAllCV <-
    glmnet::cv.glmnet(FeatAll, OutAll, lambda=costs,
                      family="binomial", alpha=0, nfolds=10)
  probFuncRvFeat <-
    getFuncRvFeat(FeatAll, logisticAllCV$glmnet.fit,
                  logisticAllCV$lambda.min)
  temp_val <-
    predict(logisticAllCV$glmnet.fit, FeatAll,
            s=logisticAllCV$lambda.min, type="response")

  ## getFuncRv_Feat works as expected
  expect_is(probFuncRvFeat, "matrix")
  expect_equal(dim(probFuncRvFeat)[1], dim(FeatAll)[1])
  expect_equal(probFuncRvFeat, temp_val)

  posteriors <-
    getFuncRvPosteriors(OutAll, probFuncRvFeat, theta=theta_init)

  ## getFuncRvPosteriors works as expected
  expect_is(posteriors, "list")
  expect_equal(rowSums(posteriors$posterior), rep(1,dim(FeatAll)[1]))

  theta.cur <- mleTheta(OutAll, posteriors$posterior, pseudoc=50)

  ## mleTheta works as expected
  expect_equal(dim(theta_init),dim(theta.cur))
  expect_equal(colSums(theta.cur),c(1,1))

  logistic.cur <- mleBeta(FeatAll, posteriors$posterior, costs)

  ## mleBeta works as expected
  beta.cur <-
    logistic.cur$beta[,which(logistic.cur$lambda == logisticAllCV$lambda.min)]
  temp_coeff <-
    glmnet::glmnet(FeatAll, posteriors$posterior,
                   lambda=logisticAllCV$lambda.min,
                   family="binomial", alpha=0)
  expect_true(all.equal(as.vector(beta.cur),
                        as.vector(temp_coeff$beta), tol=1e-04))

  logisticModelAll <-
    integratedEM(FeatAll, OutAll, logisticAllCV$lambda.min,
                 logisticAllCV$glmnet.fit, pseudoc=50, theta_init, costs)

  ## integratedEM works as expected
  expect_is(logisticModelAll, "list")
  expect_equal(logisticModelAll$lambda, logisticAllCV$lambda.min)

  postprobRIVER <- testPosteriors(FeatAll, OutAll, logisticModelAll)

  ## testPosteriors works as expected
  expect_equal(rowSums(postprobRIVER$posterior),rep(1,dim(FeatAll)[1]))
})
