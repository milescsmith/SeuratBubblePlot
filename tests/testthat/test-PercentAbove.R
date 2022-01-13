test_that("PercentAbove",{
  expect_equal(
    object = PercentAbove(x = seq(10), threshold = 5),
    expected = 50
  )
  expect_equal(
    object = PercentAbove(x = seq(10), threshold = 2),
    expected = 80
  )
  expect_equal(
    object = PercentAbove(x = seq(10), threshold = 6),
    expected = 40
  )
})