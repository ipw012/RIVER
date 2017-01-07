library(RIVERpkg)
context("Basic functions for RIVER")

test_that("Basic functions work as expected", {
  ## simulated data
  G <- simulated_features
  E <- simulated_outliers

  ## check input data
  expect_identical(G[,"SubjectID"],E[,"SubjectID"])
  expect_identical(G[,"GeneName"],E[,"GeneName"])

  g_all <- scale(as.matrix(G[,3:ncol(G)]))  # standardization
  E_disc = as.vector(as.numeric(E[,"Outlier"])) # outlier status
  theta_init=matrix(c(.99, .01, .3, .7), nrow=2)
  costs=c(100, 10, 1, .1, .01, 1e-3, 1e-4)

  ## check experimental setup
  expect_equal(length(unique(E_disc)),2)
  expect_equal(sum(unique(E_disc)),1)
  expect_equal(colSums(theta_init),c(1,1))

  cv.all.ll = glmnet::cv.glmnet(g_all, E_disc, lambda=costs,
                                family="binomial", alpha=0, nfolds=10)
  p.FR.givenG = getFRgivenG(g_all, cv.all.ll$glmnet.fit, cv.all.ll$lambda.min)
  temp_val = predict(cv.all.ll$glmnet.fit, g_all,
                     s=cv.all.ll$lambda.min, type="response")

  ## getFRgivenG works as expected
  expect_is(p.FR.givenG, "matrix")
  expect_equal(dim(p.FR.givenG)[1], dim(g_all)[1])
  expect_equal(p.FR.givenG, temp_val)

  posteriors = getFRPosteriors(E_disc, p.FR.givenG, theta=theta_init)

  ## getFRPosteriors works as expected
  expect_is(posteriors, "list")
  expect_equal(rowSums(posteriors$posterior), rep(1,dim(g_all)[1]))

  theta.cur = mleTheta(E_disc, posteriors$posterior, pseudoc=50)

  ## mleTheta works as expected
  expect_equal(dim(theta_init),dim(theta.cur))
  expect_equal(colSums(theta.cur),c(1,1))

  logistic.cur = mleBeta(g_all, posteriors$posterior, costs)

  ## mleBeta works as expected
  beta.cur = logistic.cur$beta[,which(logistic.cur$lambda == cv.all.ll$lambda.min)]
  temp_coeff = glmnet::glmnet(g_all, posteriors$posterior,
                              lambda=cv.all.ll$lambda.min, family="binomial", alpha = 0)
  expect_true(all.equal(as.vector(beta.cur), as.vector(temp_coeff$beta), tol = 1e-04))

  em.all.res <- integratedEM(g_all, E_disc, cv.all.ll$lambda.min,
                             cv.all.ll$glmnet.fit, pseudoc=50, theta_init, costs)

  ## integratedEM works as expected
  expect_is(em.all.res, "list")
  expect_equal(em.all.res$lambda, cv.all.ll$lambda.min)

  train.post = testPosteriors(g_all, E_disc, em.all.res)

  ## testPosteriors works as expected
  expect_equal(rowSums(train.post$posterior),rep(1,dim(g_all)[1]))
})
