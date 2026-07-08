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
- read_electrodePositionsOnAnnotation.m - means that we don't need getpointelectrogramname.m
- getCartoConnectorElectrodeNaming - a helper function that stores the bizarre
                                     naming conventions of catheters in Carto3!
- comparestructure.m allows two structures to be compared, 
                          eg userdata from two different versions of OpenEP

### Changed
- Surface models stored as structures rather than TriRep objects for compatibility with OpenEP
- Surface model getter and setter methods now provided
- read_ecgfile_v4.m
    - changed to read all the file data in raw form - see file
- read_meshfile.m sped up. Checks column headers explicitly as well.
- importcarto_mem.m
    - no longer has to check the electrode name for each point as this is
      taken from read_ecgfile_v4
- ability to read latest Carto data (new format without redundancy)
                
### Fixed
- Fixed issue with multiple colour bars
- Fixed issue #70 for setting tolerance when identifying in/out points in the mesh
- Geodesic distance calculator now configured for both Ubuntu and MacOSX platforms
- Tidied importcarto_mem default directories
- Sped up xml reading
- Sped up read_meshfile
- Fixed bug in getEarliestActivationSite/getLatestActivationSite when handling datasets written by EP Workbench