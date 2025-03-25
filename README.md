# EMG Analysis Ecosystem for ICMS Experiments

## Overview

This comprehensive software ecosystem provides automated analysis of electromyography (EMG) data collected during intracortical microstimulation (ICMS) experiments. The system transforms raw EMG recordings into sophisticated visualizations of muscle representation in the sensorimotor cortex through a structured, multi-level analysis pipeline.

## Key Features

- **Four-Level Analysis Architecture**: Progressively transforms data from raw signals to advanced pattern analysis
- **Automated Trial Selection**: Distinguishes valid EMG responses from artifacts using time-frequency clustering
- **Spatial Visualization**: Maps muscle responses onto cortical topography using Voronoi-based representations
- **Cross-muscle Analysis**: Reveals coordination patterns and functional modules using dimensionality reduction
- **Cross-case Comparison**: Identifies conserved organizational principles across multiple experimental subjects

## Repository Structure

- **Level_1/**: Raw data extraction, visualization, and expert assessment
  - `MAIN_TDT_Preanalysis_V11_Bat_Davis.m`: Raw data extraction and validation
  - `RawsignalPlotterV1.m`: EMG visualization engine
  - `Visual-Verdict.xlsx`: Framework for expert assessment

- **Level_2/**: Data integration, quality control, and spatial mapping
  - `Meta_Data_generator_Classification_Project_V6.m`: Integrated data consolidation
  - `Automatic_Outlier_Detector_V3.m`: Automated signal quality verification
  - `Meta_Data_plotter_flagger_V9_App.m`: Interactive signal quality assessment
  - `Voronoi_Generator/`: Scripts for generating Voronoi maps

- **Level_3/**: Timing analysis, multi-muscle visualization, and waveform topography
  - `Meta_Data_1V1_plotterV6.m`: Cross-muscle latency analyzer
  - `Stacked_Channel_Plotter_V3.m`: Multi-muscle coordination visualizer
  - `Voronoi_Waveform_plotter_V4.m`: Spatial waveform representation

- **Level_4/**: Advanced pattern recognition and cross-case analysis
  - `RDV22_TopoGraphic_Clustering_V9.m`: Hierarchical DTW-based clustering
  - `RDV22_Cross_muscle_HeatmapV1.m`: Cross-muscle pattern analysis
  - `RDV22_CrossCASE_HeatmapV3.m`: Cross-subject pattern analysis

- **Utils/**: Utility functions and classes
  - `MiladFuncV1.m`: Core utility class with signal processing and visualization functions
  - `Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform_CrossCase.m`: Voronoi processing utilities

## Requirements

- MATLAB R2019b or later
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Wavelet Toolbox (for certain Level 4 analyses)

## Installation

1. Clone this repository to your local machine
2. Add all folders and subfolders to your MATLAB path
3. See `GETTING_STARTED.md` for tutorials and example workflows

## Performance

- Reduces analysis time from 40-60 hours to 10-15 minutes per experiment
- Achieves 90-95% usable data yield compared to 55-65% with conventional methods
- Shows high correlation with expert analysis (r=0.88-0.96) for extracted features
- Successfully implemented in both rat and bat experiments, demonstrating cross-species applicability

## Citation

If you use this software in your research, please cite:
Hafezi, M., Liggins, J.A., Grochowski, L.A., et al. (2024). Topography of muscle coordination in sensorimotor cortex. [Journal details pending]
