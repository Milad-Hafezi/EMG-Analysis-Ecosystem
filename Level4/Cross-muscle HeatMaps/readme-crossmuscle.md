# Cross-muscle Heatmap (RDV22_Cross_muscle_HeatmapV1.m)

## Description
The Cross-muscle Heatmap program identifies shared activation patterns across different muscles using dimensionality reduction techniques. This MATLAB program implements non-negative matrix factorization (NMF) to discover common waveform components that may represent coordinated control signals from the motor cortex, and maps their distribution across the cortical surface using intensity-based visualization.

## Inputs
- MetaData file (.mat) from earlier processing stages
- TopoClustering output structure containing refined signatures and cluster information
- Voronoi cell definitions (processed by Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform)
- MiladFuncV1 class file containing analysis functions

## Outputs
- Comprehensive PDF documentation showing:
  - Spatial distribution maps for each component (intensity-based visualization)
  - Waveform plots showing the temporal pattern of each component
  - Combined visualizations correlating spatial patterns with waveform characteristics
- All outputs are saved in the "CrossMuscles" directory

## Key Analysis Features
1. **Non-negative Matrix Factorization (NMF)**:
   - Automatic determination of optimal number of components
   - Multiplicative update algorithm for stable factorization
   - Component waveform extraction and weight calculation
   - Reconstruction error minimization

2. **Component Attribution Analysis**:
   - Muscle-specific component scoring
   - Site-specific component mapping
   - Normalized scoring for consistent comparison
   - Handling of complex site-to-cluster relationships

3. **Spatial Visualization**:
   - Intensity-based color mapping (green gradient)
   - Component-specific cortical distribution visualization
   - Consistent visualization format across components
   - Voronoi cell-based representation

4. **Signal Processing**:
   - Optional preprocessing (mean subtraction, blanking)
   - Signal normalization
   - Time window segmentation
   - Active signal selection

## Usage
1. Run the Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform script first
2. Run the RDV22_Cross_muscle_HeatmapV1.m script
3. Select the MetaData file when prompted
4. Select the TopoClustering output file when prompted
5. The program will perform NMF analysis and generate outputs in the "CrossMuscles" directory

## Processing Pipeline
1. Load and integrate data from MetaData and TopoClustering components
2. Extract and preprocess refined signature waveforms
3. Segment the data into relevant time windows
4. Separate active and inactive signals
5. Determine optimal factorization rank
6. Perform NMF to extract components and weights
7. Calculate component scores for muscles and sites
8. Generate spatial visualizations and component waveform plots
9. Compile comprehensive PDF documentation

## Dependencies
- MATLAB (R2018b or later recommended)
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- MiladFuncV1 class file
- Voronoi cell definitions from earlier processing

## Notes
- NMF is particularly suitable for EMG analysis as it naturally enforces non-negativity
- The optimal number of components is determined automatically but limited to 20
- Intensity-based visualization provides intuitive interpretation of component strength
- This analysis reveals potential shared control signals that may coordinate multiple muscles
- The spatial mapping identifies cortical "hotspots" for each control component