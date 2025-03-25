# Topographic Clustering (RDV22_TopoGraphic_Clustering_V9.m)

## Description
The Topographic Clustering program identifies recurring waveform patterns across the cortical surface using advanced signal processing and machine learning techniques. This MATLAB program implements dynamic time warping (DTW) based hierarchical clustering to discover muscle activation patterns that may represent functional modules within the motor cortex, revealing how the brain organizes muscle coordination.

## Inputs
- MetaData file (.mat) from earlier processing stages
- Voronoi cell definitions (processed by Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform)
- MiladFuncV1 class file containing analysis functions

## Outputs
- Comprehensive Topo_Muscle_Cluster structure (.mat) containing:
  - Multi-level cluster assignments (raw, refined, latency-adjusted)
  - Activation latency data for each stimulation site
  - Cluster membership information and representative waveforms
  - Cross-channel signature patterns and their cortical distributions
  - Complete metadata for reproducibility
- PDF documentation showing:
  - Hierarchical cluster trees
  - Waveform plots for each identified cluster
  - Spatial mapping of clusters on the cortical surface
  - Comparisons between original and refined clusters

## Key Analysis Features
1. **Hierarchical DTW-based Clustering**:
   - Recursive binary splitting with configurable depth and minimum cluster size
   - DTW distance metrics that accommodate temporal variations
   - Tree visualization of cluster relationships

2. **Multi-stage Refinement Process**:
   - Small cluster reassignment (consolidates minor patterns)
   - Latency detection using multiple methods (double-threshold, adaptive, wavelet)
   - Latency-based refinement (separates timing-related patterns)
   - Normalized signal comparison (focuses on shape rather than amplitude)

3. **Spatial Mapping and Visualization**:
   - Color-coded Voronoi cells showing cluster distributions
   - Consistent color coding across visualizations
   - Before/after comparisons showing refinement benefits

4. **Cross-Muscle Integration**:
   - Processing of all muscle channels (up to 16)
   - Alignment of patterns across channels
   - Identification of coordinated patterns across muscles

## Usage
1. Run the Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform script first
2. Run the RDV22_TopoGraphic_Clustering_V9.m script
3. Select the MetaData file when prompted
4. The program will process all channels and generate outputs in the "MiladClustering" and "Topoclustering_output" directories

## Configuration Options
- **Plotters**: Controls visualization generation (1=on, 0=off)
- **maxDepth**: Maximum depth for hierarchical clustering (default: 5)
- **minClusterSize**: Minimum cluster size to allow further splitting (default: 3)
- **maxWarpingWindow**: Maximum warping window for DTW (default: 150)
- **MemberLimit**: Threshold for small cluster reassignment (default: 2)

## Processing Pipeline
1. Load and preprocess EMG data
2. Separate active and inactive sites
3. Perform hierarchical clustering on active sites
4. Reassign small clusters to larger ones
5. Detect activation latencies
6. Refine clusters based on latency information
7. Map clusters to cortical surface
8. Generate comprehensive documentation and data structures

## Dependencies
- MATLAB (R2018b or later recommended)
- Signal Processing Toolbox
- Wavelet Toolbox
- Statistics and Machine Learning Toolbox
- Parallel Computing Toolbox (for faster DTW computation)
- MiladFuncV1 class file
- Voronoi cell definitions from earlier processing

## Notes
- The clustering process may take significant time, especially for all 16 channels
- DTW-based distance calculation is computationally intensive
- For optimal results, use pre-filtered data with outliers removed
- The PDF documentation provides comprehensive visualization of the clustering results
- The Topo_Muscle_Cluster structure contains all data needed for further analysis