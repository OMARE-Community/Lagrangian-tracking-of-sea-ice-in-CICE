# README of Lagrangian Tracking in CICE (version 5)

This is a short introduction of the Lagrangian tracking in CICE (version 5) and the Community Earth System Model (CESM, version 2). The FORTRAN source code of the Lagrangian tracking module is provided, with the usage guide in the latter part of this README file. The source code of CESM is available at: [Github](https://github.com/ESCOMP/CESM). The scientific publication related to this version of Lagrangian tracking is Ning, Xu et al. (GMD-Discussions, 2024), available at: XXX. 

Prerequisite: CESM (version 2) build and running environment, including the case to be run with the Lagrangian tracking.

Usage:
- (1) Configure the case for Lagrangian tracking
- (2) Copy the provided source files in the directory of: SourceMods/src.cice
- (3) Modify the source codes if different Lagrangian tracking parameters [see also Ning, Xu et al., (GMD-Discussions, 2024)]
- (4) Compile (or re-compile) the case 
- (5) Run the case and check the case’s “run” directory for Lagrangian tracking results (a single log file for each parallel process of CICE)
