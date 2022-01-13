
test_that("make_color_scale works on viridis",{
  expected_viridis <-
    c(
      "#440154FF","#482878FF","#3E4A89FF",
      "#31688EFF","#26828EFF","#1F9E89FF",
      "#35B779FF","#6DCD59FF","#B4DE2CFF","#FDE725FF"
    )
  expect_setequal(
    object =
      make_color_scale(
        palette = "viridis",
        gradations = 10
        ),
    expected = expected_viridis
    )
})

test_that("make_color_scale works on set1",{
  expected_set1 <-
    c(
      "#E41A1C","#449B75","#AC5782",
      "#FFE528","#C66764","#999999"
    )
  
  expect_setequal(
    object =
      make_color_scale(
        palette = "Set1",
        gradations = 6
      ),
    expected = expected_set1
  )
  
})

test_that("make_color_scale works on two colors",{
  expected_red_blue <-
    c(
      "#FF0000","#CC0033","#990066",
      "#650099","#3200CC","#0000FF"
    )
  
  expect_setequal(
    object =
      make_color_scale(
        palette = c("red","blue"),
        gradations = 6
      ),
    expected = expected_red_blue
  )
  
})

test_that("make_color_scale fails on unknown palette",{
  expect_error(
    object =
      make_color_scale(
        palette = "reds",
        gradations = 6
        )
    
  )
})