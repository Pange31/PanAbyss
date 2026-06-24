# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [1.4.0]
### Added
- Checkbox to display nodes names
- Add the fasta downloading in sequence page
- Add the weighted option for local tree computation
- Add exon and transcripts information when clicking on a node
- Add p-value and score in the shared discovery regions page
- Add nodes size scale on home page
- Option for launch.sh script to generate csv import files from gfa files without GUI
- Add filter on find region discovery page
- Add a parameter to block or limit shared region discovery (for server purpose)
- Add a script for systemd configuration.
- Add a cache system.
- Add selection / coloration of genes label.
- Add the possibility to create import files or database on a cluster with apptainer

### Changed
- Gradient color and node size reference in the legend
- Round rectangle used for nodes instead of circle (configurable)
- Graph compression: option to compress only nodes whose flow exceeds a theshold
- Max nodes to get from database and max nodes to visualize are in conf file and can be changed locally on the about page.

### Fixed
- GFF annotation loading correction: the last annotation of each chromosome was associated to the next chromosome.
- Fix on the zoom functionality: when the reference genome was not present in the selection the zoom wasn't work.

### Removed


## [1.3.0]
### Added
- Store result of shared region discovery and phylogenetics in a sqlite database
- Add a limitation for shared regions discovery in conf file
- Progression bar in shared region discovery
- Add the version of PanAbyss in the about page

### Changed
- Performance improve in shared region discovery
- Change in method for retrieving annotations in shared region discovery
- Parameter the size of region where annotations are searched in the shared regions discovery page


### Fixed
- For server mode: update the multi-user functionnality

### Removed


## [1.2.0]
### Added
- Add warning message if annotation file has no common chromosome with graph
- Add legend in the visualization panel
- Add a graph for shared nodes discovery page

### Changed
- Performance update to load big annotation files
- Time consuming processes have been modified to be asynchron : shared regions discovery, global tree computation, database and annotations loading

### Fixed
-

### Removed
- 


## [1.1.0]
### Added
- Add full name in annotations
- Add a slider to set colored paths size
- Add a preset layout to allow to display more nodes if required
- Add graph compression: compact simple nodes with only one predecessor and one successor

### Changed
- Update the style of sequence output in the Shared regions discovery
- Displayed Nodes size
- Gene search has been replaced by a feature search: it is possible to select the feature to search and its value
- Add a column to get sequence into shared region discovery page
- Performance update when displaying graph of large diversity

### Fixed
-

### Removed
- 

---