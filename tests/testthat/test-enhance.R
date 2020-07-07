test_that("override of escalation::phase1_sim is nearly identical to original", {
  override_body <- body(precautionary:::phase1_sim)
  original_body <- body(escalation:::phase1_sim)
  override_lines <- deparse(override_body)
  original_lines <- deparse(original_body)
  # Of course, tests such as the following are inherently fragile.
  # But they make the point that differences are small, and it is
  # nearly 'inconceivable' that they won't fail (thus alerting me)
  # upon any changes to the original code in package 'escalation'.
  expect_equal(length(override_lines) - 2, length(original_lines))
  override_lines_omit_u_i <- grep("u_i", override_lines, invert = TRUE, value = TRUE)
  original_lines_omit_rbinom <- grep("rbinom", original_lines, invert = TRUE, value = TRUE)
  expect_equal(override_lines_omit_u_i, original_lines_omit_rbinom[-44])
})
