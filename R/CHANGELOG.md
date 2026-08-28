# Changelog

## [0.1.0] - 2026/08/28

### Added
- This changelog
- `concat` method to merge two `OlinkAssay` objects and optionally perform batch 
 correction when merging

### Changed
- Revised the default `OlinkAssay` constructor (which essentially just wraps the SummarizedExperiment constructor)
 and renamed the previous one to `OlinkAssayFromNPX`
- moved `readFromDisk` function to the `io.R` submodule

[0.1.0]: https://github.com/milescsmith/scorphan/releases/tag/v0.0.1