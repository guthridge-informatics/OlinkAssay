# Changelog

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

[0.2.0]: https://github.com/milescsmith/scorphan/releases/compare/0.1.0..0.2.0
[0.1.0]: https://github.com/milescsmith/scorphan/releases/tag/v0.1.0