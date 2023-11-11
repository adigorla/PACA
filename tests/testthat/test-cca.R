test_that("CCA works", {
  load(test_path("testdata", "cca-sims.Rdata"))

  # test 1
  cca.res.1 <- cca(cc_test1$X, cc_test1$Y, info=2)
  expect_equal(cca.res.1$corr[1:length(cc_test1$ccs)], cc_test1$ccs, tolerance = 2e-2)


  # test 2
  cca.res.2 <- cca(cc_test2$X, cc_test2$Y)
  expect_equal(cca.res.2$corr[1:length(cc_test2$ccs)], cc_test2$ccs, tolerance = 5e-2)

  # test 3
  cca.res.3 <- cca(cc_test3$X, cc_test3$Y, info=0)
  expect_equal(cca.res.3$corr[1:length(cc_test3$ccs)], cc_test3$ccs, tolerance = 9e-2)

})

