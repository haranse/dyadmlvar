context("Check lasso results")
library(dyadmlvar)

ntwrk <- read_network(fit1)

full_data <- merge(ntwrk$all, sat, by = "ID", all=TRUE, sort=TRUE)
inter_vars <- get_names(ntwrk, part = dyadmlvar::INTER, time = dyadmlvar::ALL_NETWORK)
test_res <- lasso(full_data, inter_vars, "W_csi_resid")


test_that("Correct class output", {
  expect_is(ntwrk,"dyadNetwork")
})



test_that("Correct number of names", {
  expect_equal(length(inter_vars),48)
})


test_that("Correct merge", {
  expect_equal(length(full_data),129)
  expect_equal(length(full_data$W_csi_resid),80)
})



test_that("Correct class output", {
  expect_is(test_res,"list")
  expect_equal(length(names(test_res)),2)
  expect_equal(names(test_res)[2],"Rsquared")
  expect_equal(test_res$Estimates$Importance[1],"2")
})


test_that("Correct lasso calculations", {
  expect_equal(test_res$Rsquared,0.05866679)
})
