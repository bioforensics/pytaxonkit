# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.7.2] 2021-01-27

### Added
- New `filter()` function supporting new `taxonkit filter` command (#4, #8)
- New `lca()` function supporting new `taxonkit lca` command (#5)
- Support for new `prefix` related flags (#6)
- Use new `--show-lineage-ranks` to create new `FullLineageRanks` column in `lineage()` output (#8)

### Fixed
- Bug with empty inputs (#19)


## [0.6.1] 2020-10-30

### Added
- New `name()` function taking advantage of new `--show-name` and `--no-lineage` flags in `taxonkit lineage` (#2)


## [0.6] 2020-06-12

Initial release!
