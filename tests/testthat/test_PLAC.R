context("PLAC fit")

test_that("PLAC() calls the right working horse function", {
  dat = sim.ltrc(n = 50, Cmax = 5)$dat
  expect_output(PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z1 + Z2,
                      ltrc.data = dat, td.type = "none"),
                "Calling PLAC_TI()...")
  dat = sim.ltrc(n = 50, time.dep = TRUE,
                distr.A = "binomial", p.A = 0.8, Cmax = 5)$dat
  expect_output(PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
                     ltrc.data = dat, td.type = "independent",
                     td.var = "Zv", t.jump = "zeta"),
                "Calling PLAC_TD()...")
  dat = sim.ltrc(n = 50, time.dep = TRUE, Zv.depA = TRUE, Cmax = 5)$dat
  expect_output(PLAC(ltrc.formula = Surv(As, Ys, Ds) ~ Z,
                     ltrc.data = dat, td.type = "post-trunc",
                     td.var = "Zv", t.jump = "zeta"),
                "Calling PLAC_TDR()...")
})
