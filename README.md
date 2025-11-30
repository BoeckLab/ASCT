# Introduction

This repository contains the demonstration pipeline for Antimicrobial Single-Cell 
Testing (ASCT) analysis. The demonstration pipeline can be found in the 
`Nature_microbiology_DemoPipeline` directory.

## Data Availability
Because the raw experimental data are extremely large (approximately 7 TB per 
experiment, with multiple experiments conducted), it is not feasible to share the 
complete dataset via GitHub. To address this, we provide a demonstration pipeline 
that illustrates, step by step, how a single field image and a single condition 
are processed and analyzed.

The ASCT image analysis algorithms and image data are substantial in size and 
cannot be hosted directly on GitHub. For peer review, they are freely accessible 
via Switchdrive at the following link:

[Data Access via Switchdrive](https://drive.switch.ch/index.php/s/HBJWwTGwqe2AhWW)

As an alternative, the same demonstration pipeline can also be accessed via Zenodo 
using the following private preview link:

[Data Access via Zenodo]([https://doi.org/10.5281/zenodo.172327](https://doi.org/10.5281/zenodo.17232777)

To ensure long-term accessibility, we are also preparing a deposition of these 
files in a permanent repository (Zenodo), where they will be assigned a DOI once 
the dataset is formally published. The DOI will be provided here once available.

Users may download the dataset and algorithms from these links to reproduce the 
analyses described in the demonstration pipeline.

# Repository Structure

- `1_SciCore` – Code used for high-performance computing (HPC) on the SciCore cluster.  
- `2_DataAnalysis` – Scripts used for data analysis within the HPC cluster in SciCore.  
- `3_ASCT_Figs` – Code used to render figures presented in the study.  
- `Nature_microbiology_DemoPipeline` – A demonstration pipeline showing the full 
  workflow on a reduced dataset.  

# Notes

- The provided demo pipeline is intended to illustrate the analysis steps in a 
  reproducible manner using a manageable subset of the data.  
- Please refer to the file 
  **NMICROBIOL-25031020A_nt_Code_and_software_submission_checklist.pdf** for further 
  details on code and software requirements.  
