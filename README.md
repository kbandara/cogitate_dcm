## DCM analysis of the Cogitate MEG dataset

> ***Kav Bandara** University of Melbourne, 2025*


This repo contains the analysis pipeline used for dynamic causal modelling (DCM) of the MEG dataset released by the [Cogitate Consortium (2025)](https://doi.org/10.1038/s41586-025-08888-1)).

Scripts are split into two batches (`_b1`, `_b2`), corresponding to the discovery and validation data split used in the current analysis. Note also that some scripts are functions which accept a participant index and time-window index, enabling parallelisation using job arrays on a HPC cluster (details below).


## Pipeline

The analysis has four main steps:

### 1. MEG Preprocessing (`s01_prepro`)

### 2. Source Reconstruction and ROI Localisation (`s02_source_analysis`)

### 3. Dynamic Causal Modelling (`s03_dcm`)

### 4. PEB Hypothesis Testing  (`s04_peb`)

## Dependencies

This analysis pipeline was run using matlab 2024b and SPM12. 

SPM12 is freely available from the [Wellcome Centre for Human Neuroimaging](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/).


## HPC Usage

Functions in `s01_prepro`, `s03_dcm_etw`, and `s04_dcm_offset` accept integer arguments indexing the participant list and the time window. These are designed to be submitted as SLURM jobs to run the analysis in parallel for the computationally intensive steps, e.g.:

```bash
# Submit one job per participant and time window 
sbatch --array=1-480 dcm_etw.sh  # 48 participants x 10 time windows
```

A testing block at the top of each function (commented out by default) allows running a single participant/time window within the matlab editor:

```matlab
% time_window = 1; p_names_idx = 1;
```


## Data Availability

The Cogitate MEG dataset is publicly available. Please refer to the [official data release documentation](https://cogitate-consortium.github.io/cogitate-data/). 


