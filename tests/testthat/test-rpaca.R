test_that("randomizedPACA w/ fixed K works", {
  load(test_path("testdata", "paca-sims.Rdata"))
  for ( simID in 1:5 ){
    rpaca.fixedk.res <- paca_r(powerSim[[simID]]$X, powerSim[[simID]]$Y, powerSim[[simID]]$K-1, niter = 10, batch = 800)
    rpaca.fixedk.res.estcorr <- cor(as.numeric(powerSim[[simID]]$Label),
                                    rpaca.fixedk.res$x[,2])
    cat("\n\tEst: ", rpaca.fixedk.res.estcorr,
        "\n\tTarget: ", powerSim[[simID]]$powMin , "\n")

    # check is well-powered
    expect_gte(
      abs(rpaca.fixedk.res.estcorr),
      powerSim[[simID]]$powMin
    )
  }
})
