library(matrixStats)
library(DESeq2)
library(ggplot2)

df <- arrow::read_parquet(
  "tests/initial_data/FS19030710_Extended_NPX_2026-06-22.parquet"
)
md <- readxl::read_excel(
  "tests/initial_data/2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.xlsx",
  sheet = "ManifestBuilder"
)


mat <-
  df |>
  dplyr::filter(AssayType == "assay", SampleType == "SAMPLE") |>
  dplyr::select(SampleID, Assay, Count) |>
  tidyr::pivot_wider(names_from = "Assay", values_from = "Count") |>
  tibble::column_to_rownames("SampleID") |>
  as.matrix()

md_df <- md |>
  dplyr::select(
    SampleID = `Tube ID`,
    `Sample Alias`,
    Project,
    SubjectRef,
    VisitRef
  ) |>
  dplyr::filter(SampleID %in% rownames(mat)) |>
  tibble::column_to_rownames("SampleID")

dds <- DESeqDataSetFromMatrix(countData = t(mat), colData = md_df, design = ~1)
dds <- DESeq(dds)
vsc <- vst(dds)

count_tbl <- tibble::tibble(
  variance = colVars(t(counts(dds))),
  means = rowMeans(counts(dds))
)
vsc_tbl <- tibble::tibble(
  variance = colVars(t(assay(vsc))),
  means = rowMeans(assay(vsc))
)
