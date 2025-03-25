# MetaData Generator (Meta_Data_generator_Classification_Project_V6.m)

## Description
The MetaData Generator consolidates information from Level 1 analysis into a comprehensive, structured data repository for each experimental case. It integrates verdict assessments with raw EMG data, creating a unified data structure that serves as the foundation for all higher-level analyses in the ecosystem.

## Inputs
- Excel-based verdict file (.xlsx) from the Visual-Verdict component
- MATLAB files (.mat) from the Pre-analysis code (supports multi-part experiments)
- User configuration for processing parameters

## Outputs
- Individual channel metadata files in the "MetaFiles" directory (.mat) containing:
  - Guide_trains: Matrix of EMG waveforms for all trials
  - Guide_Matrix: Matrix of metadata parameters (17 columns) for each trial
  - Muscle name, channel number, and recording parameters
  
- Combined metadata file in the "MetaFiles_Combined" directory (.mat) containing:
  - Meta_Data: Cell array of structured data for all channels
  - Each cell contains a complete data structure for one muscle channel

## Processing Features
- Signal enhancement options:
  - Notch filtering for 60Hz noise removal (configurable)
  - Butterworth low-pass filtering (adjustable parameters)
  - Artifact blanking with customizable settings
- Statistical analysis:
  - Baseline period statistics (pre-stimulation)
  - Stimulation period statistics
  - Mean and standard deviation calculations
- Threshold-based validation:
  - Integration of threshold values from verdict file
  - Classification of responses as sub-threshold or supra-threshold
  - Association of verdict status with each trial

## Usage
1. Run the Meta_Data_generator_Classification_Project_V6.m script
2. Select the verdict file (.xlsx) when prompted
3. Select data file(s) from Pre-analysis (supports multi-part experiments)
4. Configure analysis parameters:
   - Channel range to process
   - Site range to include
   - Amplitude range to process
5. The program will process all data and generate output files in the MetaFiles and MetaFiles_Combined directories

## Dependencies
- MATLAB (R2018b or later recommended)
- Excel file reading capability
- MATLAB Signal Processing Toolbox (for filtering operations)

## Hierarchical Data Organization
The program processes data through multiple organizational levels:
- Muscle channels (GOD loop)
- Data file parts (Partition loop)
- Stimulation sites (god loop)
- Stimulation amplitudes (uni loop)
- Individual trials (Chopp loop)

## Notes
- This component is a critical junction in the analysis pipeline
- The Guide_Matrix uses a standardized 17-column format including site information, stimulation parameters, timing, and verdict assessments
- The Section_no variable creates a unique identifier for each trial by combining channel, site, and amplitude information
- Both individual channel files and a combined file are created to support different analysis needs
- The output of this component feeds directly into the Outlier Detector and all higher-level analyses
