test_that("PACA w/ fixed K works", {
  load(test_path("testdata", "paca-sims.Rdata"))
  for ( simID in 1:5 ){
    paca.fixedk.res <- paca(powerSim[[simID]]$X, powerSim[[simID]]$Y, k = powerSim[[simID]]$K)
    paca.fixedk.res.estcorr <- cor(as.numeric(powerSim[[simID]]$Label),
                                   paca.fixedk.res$x[,1])
    # cat("\n\tEst: ", paca.fixedk.res.estcorr,
    #     "\n\tTarget: ", powerSim[[simID]]$powMin , "\n")
    # check is well-powered
    expect_gte(
      abs(paca.fixedk.res.estcorr),
      powerSim[[simID]]$powMin
    )
  }
})


test_that("autoPACA works", {
  load(test_path("testdata", "paca-sims.Rdata"))
  for ( simID in 1:5 ){
    paca.auto.res <- paca(powerSim[[simID]]$X, powerSim[[simID]]$Y, k = NULL)

    # check estimates approx. correct number of K components
    paca.auto.res.k <- dim(paca.auto.res$U0)[2]
    # TODO: what the expect in range function called

    # check is well-powered
    paca.auto.res.estcorr <- cor(as.numeric(powerSim[[simID]]$Label),
                                 paca.auto.res$x[,1])

    # cat("\tK_est: ", paca.auto.res.k,
    #     "\n\tEst: ", paca.auto.res.estcorr,
    #     "\n\tTarget: ", powerSim[[simID]]$powMin , "\n\n")

    expect_gte(
      abs(paca.auto.res.estcorr),
      powerSim[[simID]]$powMin
    )
  }
})


test_that("PACA null works", {
  load(test_path("testdata", "paca-sims.Rdata"))
  for ( simID in 1:3 ){
    paca.null.res <- paca_null(nullSim[[simID]]$X, nullSim[[simID]]$Y, nullSim[[simID]]$K, info=1)
    cat("\n\tNull pvalue: ", paca.null.res$pval, "\n")

    # check is calibrated
    expect_gt(
      paca.null.res$pval,
      0.05
    )
  }
})
