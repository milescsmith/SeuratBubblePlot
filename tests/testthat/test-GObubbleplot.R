library(here)
library(Seurat)
library(ggplot2)

load(here("tests/testdata/pbmc3k.RData"))
pbmc3k <- pbmc3k |>
  SCTransform() |>
  RunPCA() |>
  FindNeighbors() |>
  FindClusters() |>
  RunUMAP(
    reduction="pca",
    dims=c(1,2)
  )

save_png <- function(code, width = 400, height = 400) {
  path <- tempfile(fileext = ".png")
  png(filename = path, width = width, height = height)
  on.exit(dev.off())
  code
  path
}

test_that("test bubbleplot", {
  expect_snapshot_file(
    save_png(
      GObubbleplot(pbmc3k, go_term="GO:0009615")
    )
  )
})