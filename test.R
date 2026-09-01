x <- readFromDisk(
  npx_file = here::here(
    "tests",
    "initial_data",
    "FS19030710_Extended_NPX_2026-06-22.parquet"
  ),
  metadata_file = here::here(
    "tests",
    "initial_data",
    "2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.xlsx"
  ),
  metadata_sheet = "ManifestBuilder",
  sample_column = "Tube ID",
  verbose = TRUE
)
y <- readFromDisk(
  npx_file = here::here(
    "tests",
    "initial_data",
    "FS19030751_Extended_NPX_2026-06-22.parquet"
  ),
  metadata_file = here::here(
    "tests",
    "initial_data",
    "2026-05-15 James ABC-PreEVADE_Peds Serum_Olink (C.Guthridge)_Manifest.xlsx"
  ),
  metadata_sheet = "ManifestBuilder",
  sample_column = "Tube ID",
  verbose = TRUE
)


npx <- assay(x) |>
  dplyr::mutate(LogProtExp = ExtNPX + log2(1e5)) %>% # transform into LogProExp
  dplyr::mutate(LogProtExp_Raw = ExtNPX + log2(1e5)) # transform into LogProExp

ht_scaled_npx_sample <- npx %>%
  dplyr::filter(SampleType == "SAMPLE") %>%
  dplyr::left_join(., ht_nc_vals, by = c("Assay", "OlinkID")) %>%
  dplyr::left_join(., ht_pc_vals, by = c("Assay", "OlinkID")) %>%
  dplyr::mutate(
    sample_level_qc = dplyr::case_when(
      LogProtExp < median_nc ~ "Below LLOD",
      LogProtExp < iqr_nc ~ "Below LLOQ",
      T ~ "Pass"
    )
  ) %>%
  dplyr::mutate(
    LogProtExp = dplyr::case_when(
      LogProtExp < median_nc ~ 0,
      LogProtExp < iqr_nc ~ iqr_nc,
      T ~ LogProtExp
    )
  )


old_ht_nc_vals <- npxData %>%
  dplyr::filter(SampleType == "NEGATIVE_CONTROL") %>%
  dplyr::group_by(Assay, OlinkID) %>%
  dplyr::summarise(
    median_nc = median(na.omit(LogProtExp)),
    iqr_nc = as.numeric(quantile(na.omit(LogProtExp), 0.75)),
    .groups = 'drop'
  )

# Calculating Plate Control coefficient of variance
old_ht_pc_vals <- npxData %>%
  dplyr::filter(SampleType == "PLATE_CONTROL") %>%
  dplyr::group_by(Assay, OlinkID) %>%
  dplyr::summarise(
    pc_cv = 100 * sd(LogProtExp) / mean(LogProtExp),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(
    high_var_assay = dplyr::case_when(
      pc_cv > 20 ~ "High Variance",
      T ~ "Pass"
    )
  ) %>%
  dplyr::select(-pc_cv)

old_ht_scaled_npx_sample <-
  npxData |>
  dplyr::filter(SampleType == "SAMPLE") |>
  dplyr::left_join(ht_nc_vals, by = c("Assay", "OlinkID")) |>
  dplyr::left_join(ht_pc_vals, by = c("Assay", "OlinkID")) |>
  dplyr::mutate(
    sample_level_qc = dplyr::case_when(
      LogProtExp < median_nc ~ "Below LLOD",
      LogProtExp < iqr_nc ~ "Below LLOQ",
      T ~ "Pass"
    )
  ) |>
  dplyr::mutate(
    LogProtExp = dplyr::case_when(
      LogProtExp < median_nc ~ 0,
      LogProtExp < iqr_nc ~ iqr_nc,
      T ~ LogProtExp
    )
  )
