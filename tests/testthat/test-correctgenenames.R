library(here)
load(here("tests/testdata/processed_pbmc3k.RData"))

test_that("correctGeneNames returns HUGO names",{
  load(here("tests/testdata/corrected_gene_names.RData"))
  
  expect_setequal(
    object = rownames(x = correctGeneNames(object = pbmc3k)),
    expected = corrected_gene_names
    )
  })

test_that("correctGeneNames returns ONLY genes with HUGO name",{
  load(here("tests/testdata/filtered_corrected_gene_names.RData"))
  
  expect_setequal(
    object = rownames(
      x = correctGeneNames(
        object = pbmc3k, 
        eliminate_unknown = TRUE
        )
      ),
    expected = filtered_corrected_gene_names
    )
})

