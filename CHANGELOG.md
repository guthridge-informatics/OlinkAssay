# Changelog

## [0.3.0] - 2026/09/03

### Added
- now storing the `SampleQC` and `AssayQC` in colData and rowData, respectively as nested lists.
- added `retrieveAssayQC` and `retrieveSampleQC` methods to appropriately unnest their target data.
- `OlinkAssayFromNPX` now takes an "extra_metadata" parameter that allows for external data (i.e. lot numbers) to be added to the `metadata(obj)` list
- Stole `colData` and `rowData` replacement methods from `{SummarizedExperiment}`
- Added a `QCPlot()` function to plot sample and assay QC data.

### Changed
- Sanitizing SampleIDs on import, meaning *NO* samples IDs starting with a number and nothing with a period or dash in the name; leading numbers are prefixed with an "X" and internal symbols are converted to an underscore. Anything being converted to a S4Vector::DataFrame gets mangled if you have these in your sample IDs, leading to disconnects between data and metadata.
- Storing assays as matrices
- Added storing the "AssayQC" and "AssayQCWarn" metrics to the `rowData` slot; since each "PlateID" and the metrics are actually generated for each subplate, these are stored as a `{PlateID}:{result}` list.
- Partially revised the `performance_report_template.qmd` to accommodate using an OlinkAssay object. Still a work in progress.
- renamed `medianCorrection` to `batchCorrection`

### Fixed
- merging two OlinkAssay objects preserves row order so we don't end up with a weird assay name:OlinkID disconnect.
- overhauled how `rowData` was merged between two objects; created a `.combineRowData` function to take care of it

## [0.2.0] - 2026/08/21

### Added
- documentation
- Add two new slots to `OlinkAssay` - `negativeControls` and `plateControls`, 
 `S4Vectors::DataFrame` that store the negative and plate control data, respectively.
- Add:
 - `scaleNPX` to transform the ExtNPX data to a non-negative matrix
 - `level2Process` to handle the data standard level2 QC.

### Changed
- removed `as_SummarizedExperiment()` since it was duplicating `OlinkAssayFromNPX()`
- revised OlinkAssay creation
 - removed the `.OlinkAssay` function and now more closely mimic the `SummarizedAssay` constructor
 - began adding checks to constructor

### Fixed
- Moved CHANGELOG.md to its proper location
- extra metadata is now actually saved to the OlinkAssay object on import from disk

## [0.1.0] - 2026/08/28

### Added
- This changelog
- `concat` method to merge two `OlinkAssay` objects and optionally perform batch 
 correction when merging

### Changed
- Revised the default `OlinkAssay` constructor (which essentially just wraps the SummarizedExperiment constructor)
 and renamed the previous one to `OlinkAssayFromNPX`
- moved `readFromDisk` function to the `io.R` submodule

[0.3.0]: https://github.com/milescsmith/scorphan/releases/compare/0.2.0..0.3.0
[0.2.0]: https://github.com/milescsmith/scorphan/releases/compare/0.1.0..0.2.0
[0.1.0]: https://github.com/milescsmith/scorphan/releases/tag/v0.1.0