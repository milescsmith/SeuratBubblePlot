## code to prepare `DATASET` dataset goes here

usethis::use_data(pbmc3k, overwrite = TRUE)

library(Seurat)
load(here::here("tests/testdata/pbmc3k.RData"))
pbmc3k <- pbmc3k |> 
  SCTransform() |>
  RunPCA() |> 
  FindNeighbors() |> 
  FindClusters() |> 
  RunUMAP(
    reduction="pca",
    dims=c(1,2)
  )