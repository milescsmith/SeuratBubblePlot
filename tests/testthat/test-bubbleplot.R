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
      bubbleplot(
        object = pbmc3k, 
        features_plot=c("IFIT1","IFITM3","ISG15","IFNAR1")
        )
    )
  )
})
