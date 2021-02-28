# Phase-Shift-ST
2D+1 ST Testing and Study

This is the code underlying Wright et al (2021), a manuscript submitted to Atmospheric Measurement Techniques which applies the "2D+1" S-Transform to real and synthetic AIRS data to measure atmospheric gravity waves.


For the code to work, you will need to:
1. add the directory "functions/" to your Matlab path
2. add the top-level path for the full set of files to "functions/LocalDataDir.m". 

After this the code *should* hopefully just work (it does in my tests). Note that sections plotting topography and surface maps have been commented
out - this is because the filesizes of these datasets are too large for this medium. Data used for his purpose in the paper figures were obtained
from NaturalEarth and the USGS.


Figure 1 is generated using "05FlatvsStretch/plot_waves_3d.m"

Figure 2 is generated using "02Conceptual Diagrams/plot_conceptual_bits.m". 

Figure 3 is generated using "03SyntheticWaveTests/01AveragingVolumeDiagram/diagram.m"

Figure 4 is generated using first "03SyntheticWaveTests/02SystematicTests/generate_data_v2.m" and then "03SyntheticWaveTests/02SystematicTests/plot_comparisons.m"

Figure 5 is generated using "04RealWaveTests/01WaveImages/plot_waves_3d.m"

Figures 6-8 are generated using the files "results_example1.m", "results_example2.m" and "results_example3.m" in "04RealWaveTests/02Analyses"


The directory "01GettingItWorking" just contains files I used to test concepts, which are incomplete and are likely to be missing required functions. 
Ignore these, they are just here for completeness.


For all figures, manual modification to the layout, including adding contextual lines and rearranging panels into figures, was carried out 
manually in image-editing software. So: the above files will not produce final versions as seen in the paper - the contents should be fully 
consistent however, and demonstrate that the figures accurately represent the input data and our analysis.


The data/ directory includes the AIRS granules used as tests - these were produced by Lars Hoffman as described in Hoffman and Alexander (JGR, 2009).
 
