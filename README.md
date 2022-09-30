# Multishot-EPI
We proposed an image reconstruction method based on Hankel structured low-rank matrix completion and an interleaved 3D CAIPI sampling pattern, which can improve the robustness to physiological fluctuations of 3D multi-shot EPI acquisitions for fMRI at 7T. 

run_multishot_SLR.m is the main function, which performes the SENSE reconstruction for the standard blipped-CAIPI data, and the proposed structured low-rank reconstruction for the seg-CAIPI data. The optimization problem of the proposed reconstruction is solved by ADMM algorithm.  

A reference dataset can be found here https://doi.org/10.5281/zenodo.7128625

The simulation data in the paper is also provided
