# suski-2022-tracking
Cell tracking an analysis code used for Suski et al. 2022, Nature.

Image analysis pipeline has primarly been  used on MCF-10A cells, and has been used on HeLa and RPE-1 cells. This pipeline has the following functionaliy:
1) Tracking cells from time-lapse microscopy and quantification of fluorescent reporters (nuclear/cytoplasmic signals, foci analysis)
2) Quantitative image-based cytometry (QIBC) analysis of fixed-cell microscopy, with capabilities to handle multi-round imaging of the same sample.
3) Retrospective Time-lapse Synchronized QIBC (RT-QIBC) analysis of fixed-cell microscopy, using prior information from live-cell tracking measurements of matched cells. 

See Ratnayeke et al. 2021, BioRxiv (https://doi.org/10.1101/2021.07.06.451311) for methodological details.  Pipeline was written based on code from Cappell et al., 2016, Cell 166, 167-180. June 30, 2016. 
## Installation
Add *cell-cycle-tracking/Tracking_dependencies/* to MATLAB path. 
## Usage
In *Tracking_main/*, main functions for analysis are found for each microscope. 
