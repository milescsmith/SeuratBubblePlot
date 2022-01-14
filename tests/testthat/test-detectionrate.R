load("tests/testdata/detection_rate_expected.RData")

test_that("DetectionRate",{
  expect_equal(
    object = DetectionRate(pbmc3k, c("IFIT1", "IFNAR1","IFITM3"), assay="RNA"),
    expected = isg_rna_detect_rate
    )
  expect_equal(
    object = DetectionRate(pbmc3k, c("IFIT1", "IFNAR1","IFITM3"), assay="SCT"),
    expected = isg_sct_detect_rate
  )
})
