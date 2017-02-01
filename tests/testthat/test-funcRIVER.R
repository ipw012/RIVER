library(RIVERpkg)
context("Basic functions for RIVER")

test_that("Basic functions work as expected", {
  ## getData
  dataInput = getData(filename=system.file("extdata", "simulation_RIVER.gz",
                                           package = "RIVERpkg"), ZscoreThrd=1.5)
  GAll = t(exprs(dataInput)) # all genomic features
  EAll = as.numeric(unlist(dataInput$Outlier))-1 # all outlier status
  expect_equal(class(GAll),"matrix")
  expect_equal(class(EAll),"numeric")
  expect_equal(length(unique(EAll)),2) # either outlier or nonoutlier

  Gtrng = t(exprs(dataInput[,is.na(dataInput$N2pair)])) # training G data
  Etrng = as.numeric(unlist(dataInput$Outlier[is.na(dataInput$N2pair)]))-1 # training E data
  expect_equal(class(Gtrng),"matrix")
  expect_equal(class(Etrng),"numeric")

  Gtest = t(cbind(exprs(dataInput[,!is.na(dataInput$N2pair)])
                  [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                  exprs(dataInput[,!is.na(dataInput$N2pair)])
                  [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)])) # test G data
  # test E data (indiv1 and then indiv2 from N2 pairs)
  Etest1 = as.numeric(unlist(c(dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                               dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1
  # test E data (indiv2 and then indiv1 from N2 pairs)
  Etest2 = as.numeric(unlist(c(dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)],
                               dataInput$Outlier[!is.na(dataInput$N2pair)]
                               [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)])))-1
  expect_equal(class(Gtest),"matrix")
  expect_identical(rownames(Gtest),c(sampleNames(dataInput[,!is.na(dataInput$N2pair)])
                         [seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                         sampleNames(dataInput[,!is.na(dataInput$N2pair)])
                         [seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))
  expect_equal(length(unique(Etest1)),2) # either outlier or nonoutlier
  expect_equal(length(unique(Etest2)),2) # either outlier or nonoutlier

  # N2 pairs should have same genomic features
  expect_equivalent(sum(diag(cor(exprs(dataInput[,!is.na(dataInput$N2pair)])
                                 [,seq(from=1,to=sum(!is.na(dataInput$N2pair)),by=2)],
                                 exprs(dataInput[,!is.na(dataInput$N2pair)])
                                 [,seq(from=2,to=sum(!is.na(dataInput$N2pair)),by=2)]))),
                    sum(!is.na(dataInput$N2pair))/2)

  theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
  costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)

  ## check experimental setup
  expect_equal(sum(unique(EAll)),1)
  expect_equal(colSums(theta_init),c(1,1))

  cv.all.ll = glmnet::cv.glmnet(GAll, EAll, lambda=costs,
                                family="binomial", alpha=0, nfolds=10)
  p.FR.givenG = getFRgivenG(GAll, cv.all.ll$glmnet.fit, cv.all.ll$lambda.min)
  temp_val = predict(cv.all.ll$glmnet.fit, GAll,
                     s=cv.all.ll$lambda.min, type="response")

  ## getFRgivenG works as expected
  expect_is(p.FR.givenG, "matrix")
  expect_equal(dim(p.FR.givenG)[1], dim(GAll)[1])
  expect_equal(p.FR.givenG, temp_val)

  posteriors = getFRPosteriors(EAll, p.FR.givenG, theta=theta_init)

  ## getFRPosteriors works as expected
  expect_is(posteriors, "list")
  expect_equal(rowSums(posteriors$posterior), rep(1,dim(GAll)[1]))

  theta.cur = mleTheta(EAll, posteriors$posterior, pseudoc=50)

  ## mleTheta works as expected
  expect_equal(dim(theta_init),dim(theta.cur))
  expect_equal(colSums(theta.cur),c(1,1))

  logistic.cur = mleBeta(GAll, posteriors$posterior, costs)

  ## mleBeta works as expected
  beta.cur = logistic.cur$beta[,which(logistic.cur$lambda == cv.all.ll$lambda.min)]
  temp_coeff = glmnet::glmnet(GAll, posteriors$posterior,
                              lambda=cv.all.ll$lambda.min, family="binomial", alpha = 0)
  expect_true(all.equal(as.vector(beta.cur), as.vector(temp_coeff$beta), tol = 1e-04))

  em.all.res <- integratedEM(GAll, EAll, cv.all.ll$lambda.min,
                             cv.all.ll$glmnet.fit, pseudoc=50, theta_init, costs)

  ## integratedEM works as expected
  expect_is(em.all.res, "list")
  expect_equal(em.all.res$lambda, cv.all.ll$lambda.min)

  train.post = testPosteriors(GAll, EAll, em.all.res)

  ## testPosteriors works as expected
  expect_equal(rowSums(train.post$posterior),rep(1,dim(GAll)[1]))
})
