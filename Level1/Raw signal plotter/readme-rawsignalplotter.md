# Raw Signal Plotter (RawsignalPlotterV1.m)

## Description
The Raw Signal Plotter generates visual reports of EMG signals from ICMS experiments, organized by stimulation site, current amplitude, and muscle channel. It transforms the structured data from the Pre-analysis code into comprehensive PDF reports featuring stacked EMG traces with statistical references, enabling researchers to efficiently assess muscle responses to cortical stimulation.

## Inputs
- MATLAB files (.mat) from the Pre-analysis code containing EMG data and metadata
- Support for multi-part experiments with separate data files
- User-defined configuration settings for visualization and processing parameters

## Outputs
- PDF reports (one for each muscle channel) containing:
  - Stacked EMG traces with appropriate vertical separation
  - Filtered mean traces overlaid for comparison
  - Stimulation timing indicators (vertical lines at 0ms and 500ms)
  - Statistical reference lines (baseline mean and standard deviation)
  - Trial identification and timing information
  - Organized by stimulation site and current amplitude

## Usage
1. Run the RawsignalPlotterV1.m script
2. Select data files when prompted (supports single or multi-part experiments)
3. Configure analysis parameters:
   - Channel range to process
   - Site range to include
   - Amplitude range to visualize
   - Processing options (filtering, blanking, display parameters)
4. PDF reports will be generated in a date-stamped folder (PDFPlots_dd-mmm-yyyy)

## Dependencies
- MATLAB (R2018b or later recommended)
- MATLAB Report Generator Toolbox (for PDF creation)
- MATLAB Signal Processing Toolbox (for filtering operations)

## Processing Features
- Signal enhancement:
  - Notch filtering for 60Hz noise removal (configurable)
  - Butterworth low-pass filtering (adjustable order and cutoff)
  - Artifact blanking with customizable parameters
  - Trimmed mean calculation for robust averaging
- Visualization optimization:
  - Adaptive vertical scaling based on signal characteristics
  - Automatic handling of signals with extreme outliers
  - Clear separation between trials for visual assessment
  - Comprehensive labeling of experimental parameters

## Notes
- The generated PDFs serve as the foundation for the Visual-Verdict component
- Supports both batch and single-channel processing modes
- Includes options for processing both blanked and unblanked EMG data
- PDF reports include timestamps, case identifiers, and experimental parameters for traceability
