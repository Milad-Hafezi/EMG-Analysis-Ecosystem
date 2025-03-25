# 1v1 Latency Analyzer (Meta_Data_1V1_plotterV6.m)

## Description
The 1v1 Latency Analyzer enables detailed investigation of timing relationships between muscle activations evoked by intracortical microstimulation. This MATLAB program uses cross-correlation techniques to quantify and visualize the temporal sequence of muscle recruitment, providing insights into how the sensorimotor cortex coordinates muscle activity.

## Inputs
- MetaData file (.mat) from the MetaData Generator or outlier detection components
- User selections for:
  - Stimulation site number
  - Stimulation amplitude
  - Muscles of interest for comparison

## Outputs
- Summary plots showing:
  - Raw and filtered averaged responses for all selected muscles
  - Stimulation period markers (onset/offset)
- Matrix comparison plots showing:
  - Pairwise timing relationships between all selected muscles
  - Optimal temporal alignment between each muscle pair
  - Quantified lag values in milliseconds
  - Available in multiple analysis modes (mean/median, filtered/unfiltered)
- Optional individual trial comparison plots when enabled

## Analysis Modes
1. **Mean-based analysis with filtering**: 
   - Averages all trials for each muscle
   - Applies Butterworth filtering (20Hz cutoff)
   - Provides noise-robust timing assessment

2. **Median-based analysis with filtering**:
   - Uses median response for outlier resistance
   - Applies Butterworth filtering (20Hz cutoff)
   - Robust when spontaneous activity may be present

3. **Mean/Median without filtering**:
   - Performs timing analysis on raw signals
   - Useful for examining frequency-dependent timing relationships

4. **Trial-by-trial analysis**:
   - Examines timing relationships in individual trials
   - Useful for assessing consistency across repetitions
   - Interactive interface with waiting between trials

## Usage
1. Run the Meta_Data_1V1_plotterV6.m script
2. When prompted, enter the site number for analysis
3. When prompted, enter the stimulation amplitude
4. When prompted, enter the channel numbers to compare (separated by spaces)
5. View the summary plots showing all selected channels
6. Examine the matrix comparison plots showing timing relationships
7. (Optional) If V1V=1, view trial-by-trial comparisons

## Configuration Options
- **Plot**: Controls individual trial plotting (0=off, 1=on)
- **V1V**: Controls trial-by-trial correlation analysis (0=off, 1=on)
- Filter parameters:
  - Summary plots: 45Hz Butterworth filter
  - Comparison plots: 20Hz Butterworth filter

## Technical Features
- Cross-correlation-based timing analysis
- Zero-phase filtering to prevent timing distortion
- Comprehensive matrix visualization framework
- User-friendly channel selection with muscle name display
- Millisecond-precision timing measurements

## Dependencies
- MATLAB (R2018b or later recommended)
- Signal Processing Toolbox (for filtering and cross-correlation)
- MetaData file from earlier processing stages

## Notes
- Timing values (lags) are reported in milliseconds
- Positive lag indicates the second signal occurs after the first
- For optimal timing analysis, use preprocessed data with outliers removed
- The matrix visualization is organized with reference muscles in rows and comparison muscles in columns