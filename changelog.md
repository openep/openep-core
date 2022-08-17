# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).
 
## [Unreleased] - yyyy-mm-dd
 
This section documents changes which have been merged into develop and will form the next release
 
### Added
- Function to anonymise a directory of OpenEP datasets
- Function to convert surface models from TriRep to structure
- Function to convert surface models from TriRep to structure for a directory of OpenEP datasets
- Parsing data from the Precision electroanatomic mapping system
- Parsing data from the EnSiteX electroanatomic mapping system

### Changed
- Surface models stored as structures rather than TriRep objects for compatibility with OpenEP
- Surface model getter and setter methods now provided

### Fixed
- Fixed issue with multiple colour bars
- Fixed issue #70 for setting tolerance when identifying in/out points in the mesh
- Geodesic distance calculator now configured for both Ubuntu and MacOSX platforms