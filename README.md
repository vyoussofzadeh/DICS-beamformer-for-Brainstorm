# DICS-beamformer-for-Brainstorm
This is a DICS beamformer implementation for the Brainstorm (BS) software package. The dynamic imaging of coherent sources (DICS) beamformer technique enables the study of cortical sources of oscillatory activation in frequency-domain (Gross et al., 2001). DICS is a linearly constrained minimum variance beamformer in the frequency domain. It estimates the covariance matrix to calculate the spatial filter using the sensor-level cross-spectral density (CSD) matrix and applies the filter to the sensor-level CSD to reconstruct the source-level CSDs of pairwise voxel activations, providing coherence measures between the source pairs.

Note that the new pipeline is available at the BS GitHub repository: https://github.com/brainstorm-tools/brainstorm3/blob/master/toolbox/process/functions/process_ft_sourceanalysis_dics.m. The following details are based on the older version of the pipeline.

This implementation has mainly focused on localizing induced activations due to task-MEG responses, e.g., an overt definition naming language experiment.

Before running, follow these steps:

- Add the FieldTrip toolbox to the Matlab path, e.g., ft_path = 'xx/fieldtrip_20190419'; addpath(ft_path); ft_defaults;
- Estimate headmodel: overlapping spheres for surface-based and MRI volume for volumetric-based source mapping;

To run the DICS-BF in BS:

1. Open Brainstorm
2. Add (preprocessed epoched) trial responses to the processing window,<br/>
<p align="center">
<img src="images/1_screenshot.png" width="500">
</p>
3. Select the DICS-BF source modeling from the process selection/Source/FieldTrip: ft_souceanalysis DICS-BF, vXX, and Run, <br/>
<p align="center">
<img src="images/2_screenshot.png" width="500">
</p>
4. Choose DICS-beamformer as the source modeling approach, and MEG (MEG-MAG, or MEG GRAD) as the sensor type, 
<p align="center">
<img src="images/9_screenshot.png" width="400">
</p>
5. Pipeline estimates time-frequency responses (sensor-space, average across all sensors), <br/>
<p align="center">
<img src="images/3_screenshot.png" width="400">
</p>
6. Select the time interval of post-vs-pre responses, e.g., [-0.3,0;0.7,1.2]
7. Select the frequency of interest, e.g., f=22Hz; A dpss smoothing window of 4Hz is applied (by default, see `vy_fft`, line 656) to estimate the cross-spectral density (CSD) matrix.
8. Results (surface map) are stored in the last trial response.<br/>
<p align="center">
<img src="images/8_screenshot.png" width="400">
</p>
9. A sample result, an auditory definition naming task, DICS-BF compared against a dynamic Statistical Parametric Maps (dSPM), broadband 0.1-28Hz, is provided below. <br/>
<p align="center">
<img src="images/7_screenshot.png" width="600">
</p>
For further inquiries, please contact vyoussofzadeh@mcw.edu.
 
# Cite
1. Gross J, Kujala J, Hamalainen M, Timmermann L, Schnitzler A, Salmelin R. Dynamic imaging of coherent sources: Studying neural interactions in the human brain. Proc Natl Acad Sci U S A. 2001;98(2):694–9.
2. Youssofzadeh, V., Stout, J., Ustine, C., Gross, W.L., Lisa, L., Humphries, C.J., Binder, J.R., Raghavan, M., 2020. Mapping language from MEG beta power modulations during auditory and visual naming, NeuroImage. Elsevier Inc. https://doi.org/10.1016/j.neuroimage.2020.117090

# Updates
- 07/08/21, A new tutorial wiki BS page was created: https://neuroimage.usc.edu/brainstorm/Tutorials/DICS
- 07/08/21, Pipeline is officially added to BS github repository: https://github.com/brainstorm-tools/brainstorm3/blob/master/toolbox/process/functions/process_ft_sourceanalysis_dics.m
- 07/01/21, a new update was made. Script was reformatted to match the standards in Brainstorm.
- 06/06/21, a new version of the pipeline [process_ft_sourceanalysis_dics.m](https://github.com/vyoussofzadeh/DICS-beamformer-for-Brainstorm/blob/master/process_ft_sourceanalysis_dics.m) was added to the repository. For convenience and compatibility, variable inputs were integrated into a BS option GUI interface. The new version should produce the same results as the old version [process_ft_sourceanalysis_DICS_BF.m](https://github.com/vyoussofzadeh/DICS-beamformer-for-Brainstorm/blob/master/Older%20version/process_ft_sourceanalysis_DICS_BF.m).
   

<!-- <img src="images/4_screenshot.png" width="500"> -->
<!-- <img src="images/5_screenshot.png" width="600"> -->
<!-- <img src="images/6_screenshot.png" width="600"> -->
