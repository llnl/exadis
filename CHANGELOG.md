# Change Log
All notable changes to ExaDiS will be documented in this file.


## Version 0.2.0 - Mar 4, 2026

### Added
- Added registry mechanisms for force and mobility types
- Added pybind submodule registry mechanism
- Added a more flexible ForceGlobal model
- Added SYCL compilation option
- Added crystal option num_bcc_plane_families to allow for 110, 112, and 123 BCC planes
- Added check for NaNs in time-integrators
- Added build option to use compute maps instead of atomic operations
- Added force map option in ForceFFT and ForceSeg classes
- Added binding to internal functions for manipulating SerialDisNet instances
- Added new coarsening mode in Remesh class


### Changed
- Changed FCC_0_FRIC mobility to have offset and scaling factors
- Refactored force binding for more direct access to backend classes
- Added check for crystallographic planes when testing for violations in collision
- Modified the collision function to better handle hinges
- Allowed for variable number of subgroups in subcycling integration


### Fixed
- Circumvented neighbor list miscompilation bug with CUDA 12.6+
- Removed case-sensitive check for Kokkos_ENABLE_CUDA option
- Added missing fence in test_kokkos
- Checking for use of exadis_plastic_strain in pyexdis driver 


## Version 0.1.4 - Mar 31, 2025

### Added
- Added FCC_0_FRIC mobility law with spatial field option
- Added cross-slip module for FCC and BCC crystals
- Added ExaDisNet utility method generate_line_config()
- Added crystal orientation option to ExaDisNet generate methods
- Extended ForceFFT to work with triclinic boxes
- Added mobility law FCC_0B
- Added optional friction stress to BCC_0B mobility law
- Added option to use glide planes for BCC crystals
- Added the OpRec module
- Added ExaDiS tools to compute dislocation fields
- Added helper function to manage resizing of large Kokkos views


### Changed
- Split the build of the python module to reduce compilation times on devices
- Refactored TopologyParallel to make the SplitDisNet class more widely available
- Refactored crystal options and added parameter enforce_glide_planes
- Optimized implementation of FCC_0 mobilities
- Refactored the implementation of the integrators
- Refactored the implementation of the ForceSegSegList/SegSegList classes
- Optimized the neighbor list classes
- Optimized various kernels and launch bounds for AMD GPUs
- Upgraded to Kokkos 4.6.00


### Fixed
- Fixed bug for resetting matrices at restart
- Fixed memory leak in ExaDisNet binding class


## Version 0.1.3 - Aug 19, 2024

### Added
- Added force cutoff model and corresponding bindings
- Extended ForceFFT to support different resolution in each direction
- Added compute_team option to ForceSeg and team kernel to ForceN2
- Added binding to the simulation stepping options
- Added binding to the simulation restart option
- Added option to automatically select the minseg param
- Added option strict to NeighborBox class
- Added binding to enable calling python-based modules from within ExaDiS
- Added benchmark tests
- Added CHANGELOG file


### Changed
- Added function register_neighbor_cutoff()
- Modified binding to accept independent FFT grid dimensions
- Switched to using the Voigt convention for stress in the python binding to be consistent with OpenDiS


### Fixed
- Checking for positive pair distances when building seg/seg lists


## Version 0.1.2 - Jun 22, 2024

### Changed
- Made modifications to allow for hybrid use of unified memory

### Fixed
- Fixed bug to use shared memory for armsets array in TopologyParallel
- Fixed bugs with the use of triclinic cells and modified Cell to use origin instead of min-max coords


## Version 0.1.1 - Jun 4, 2024

### Added
- Added use of node tags internally and in python binding
- Added tests for kokkos and cuda
- Added non-linear BCC mobility law Mobility_BCC_nl
- Added helper functions make_system() make_network_manager()
- Added remesh rule to not refine between pinned nodes

### Changed
- Modified test cases to work with updated DisNetManager
