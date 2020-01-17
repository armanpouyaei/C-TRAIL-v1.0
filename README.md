C-TRAILv1.0 is based on CMAQv5.2 (you can find all the changes in the /CCTM/scripts/ and /CCTM/src/)
========
If you are planning to run the C-TRAIL model, you need to use CMAQ's tutorials.
========

C-TRAIL/CMAQ combines current knowledge in atmospheric science and air quality modeling with multi-processor
computing techniques in an open-source framework to deliver fast, technically sound estimates of ozone,
particulates, toxics, and acid deposition.
## C-TRAIL Repository Guide
Source code and scripts are organized as follows:
* **CCTM:** code and scripts for running the 3D-CTM at the heart of C-TRAIL.
* **DOCS:** CMAQ User's Guide, developers guidance, tutorials and known issues.
* **PREP:** Data preprocessing tools for important input files like initial and boundary conditions, etc.
* **POST:** Data postprocessing tools for aggregating and evaluating CMAQ output products (e.g. Combine, Site-Compare, etc)
* **TUTORIALS:** Short tutorials provide examples of different CMAQ features.
* **UTIL:** Utilities for generating code and using CMAQ (e.g. chemical mechanism generation, IO-API, etc)

## Documentation
Release notes and Code documentation are included within this repository (they are version-controlled along with the code itself).  

[CMAQv5.2 Documentation](DOCS/README.md)   
[CMAQv5.2 Release Notes](CCTM/docs/Release_Notes/README.md)   
[CMAQv5.2 Known Issues](DOCS/Known_Issues/README.md)   

## CMAQ Test Cases
Benchmark/tutorial input and output data for each CMAQ release version are available from the CMAS Center Software Clearinghouse. From http://www.cmascenter.org, select Download -> Software -> CMAQ and choose the version to get the tutorial data.  

