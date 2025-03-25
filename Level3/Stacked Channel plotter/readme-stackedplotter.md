# Stacked-muscle Plotter (Stacked_Channel_Plotter_V3.m)

## Description
The Stacked-muscle Plotter generates comprehensive visualizations displaying simultaneous EMG activity across all recorded muscles in response to stimulation at a specific cortical site. This MATLAB program creates stacked, time-aligned displays that reveal muscle coordination patterns evoked by intracortical microstimulation, providing insight into functional muscle synergies.

## Inputs
- MetaData file (.mat) from the MetaData Generator or outlier detection components
- User selections for:
  - Stimulation site number
  - Stimulation amplitude

## Outputs
- Comprehensive stacked visualization showing:
  - All muscle responses to the selected site and amplitude
  - Time-aligned traces referenced to stimulation onset/offset
  - Clear vertical separation between muscles
  - Channel numbers and threshold values for each muscle
- Output formats:
  - PDF reports (when PDF_logic=1)
  - SVG files for publication-quality vector graphics

## Configuration Options
- **Separation_type**: Method for calculating vertical spacing ('Tec' uses signal amplitude characteristics)
- **Separation**: Ratio for vertical spacing between traces (multiplier for calculated separation)
- **Stacked_filtration**: Enables pre-averaging signal filtering (0=off, 1=on)
- **Post_filt**: Enables post-averaging signal filtering (0=off, 1=on)
- **Butter_Cut_Off**: Cutoff frequency for Butterworth filtering (default: 45Hz)
- **Order**: Filter order for Butterworth filter (default: 4)
- **Normalized_Train**: Toggles amplitude normalization (0=off, 1=on)
- **PDF_logic**: Controls PDF report generation (0=off, 1=on)

## Usage
1. Run the Stacked_Channel_Plotter_V3.m script
2. If not already loaded, select the MetaData file when prompted
3. When prompted, enter the site number for visualization
4. When prompted, enter the stimulation amplitude
5. The program will generate a stacked visualization of all muscle responses
6. Output files will be saved in the "Stacked_muscles" and "Stacked_musclesSVG" directories

## Technical Features
- Trimmed mean calculation (40% trimming) for robust averaging
- Configurable pre- and post-averaging filtering
- Intelligent vertical separation based on signal characteristics
- Comprehensive labeling with experimental parameters
- Time-aligned visualization relative to stimulation onset/offset
- Support for both normalized and raw amplitude display

## Signal Processing
- **Trimmed Mean**: Uses trimmed mean (rejecting highest and lowest 20% of values) for robust averaging
- **Butterworth Filtering**: Implements low-pass filtering with configurable cutoff and order
- **Normalization**: Optional amplitude normalization for consistent display across muscles

## Dependencies
- MATLAB (R2018b or later recommended)
- Signal Processing Toolbox (for filtering)
- Report Generator Toolbox (for PDF generation when PDF_logic=1)
- MetaData file from earlier processing stages

## Notes
- The visualization is time-aligned with vertical dashed lines marking stimulation onset (0ms) and offset (500ms)
- Channel numbers and threshold values are displayed to the right of each trace
- For publication-quality figures, SVG format is recommended
- The configuration options can be adjusted at the top of the script before running