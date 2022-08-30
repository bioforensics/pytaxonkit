# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [0.8] 2022-08-30

## Changed
- Updated doctest for `pytaxonkit.list` to include taxa that are updated in the NCBI taxonomy DB less frequently (#13, #14, #21, #24)
- Other occasional updates to compensate for changes in the NCBI taxonomy DB (#15, #16, #17, #23, #26)
- Updated `pytaxonkit.list_ranks()` to match updated `taxonkit filter --list-order` behavior (#20)
- Now using Black instead of pycodestyle to check and autoformat Python code (#27)

### Fixed
- Bug with how certain functions handle empty inputs (#19)
- Bug causing `pytaxonkit.lineage` to fail when only a single taxid was provided as input (#22)
- Bug causing `pytaxonkit.lca` to fail when all input taxids are deleted or "unfound" (#22)
- Bug causing `pytaxonkit.name` to fail on empty inputs (#26)


## [0.7.2] 2021-01-27

### Added
- New `filter()` function supporting new `taxonkit filter` command (#4, #8)
- New `lca()` function supporting new `taxonkit lca` command (#5)
- Support for new `prefix` related flags (#6)
- Use new `--show-lineage-ranks` to create new `FullLineageRanks` column in `lineage()` output (#8)


## [0.6.1] 2020-10-30

### Added
- New `name()` function taking advantage of new `--show-name` and `--no-lineage` flags in `taxonkit lineage` (#2)


## [0.6] 2020-06-12

Initial release!
