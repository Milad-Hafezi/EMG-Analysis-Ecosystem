# Voronoi Waveform Plotter (Voronoi_Waveform_plotter_V4.m)

## Description
The Voronoi Waveform Plotter creates spatial visualizations where actual EMG waveforms are directly mapped onto their corresponding stimulation sites across the cortical surface. This MATLAB program combines EMG data with Voronoi cell definitions to provide a comprehensive view of how muscle response characteristics vary topographically across the sensorimotor cortex.

## Inputs
- MetaData file (.mat) from the MetaData Generator or outlier detection components
- SVG file containing Voronoi cell definitions (processed by Milad_Davis_code)
- Configuration parameters for visualization and signal processing options

## Outputs
- Comprehensive spatial visualization (SVG format) showing:
  - Complete cortical map with Voronoi cell boundaries
  - EMG waveforms overlaid at each stimulation site
  - Highlighted stimulation periods
  - Site labels with activity markers
  - Scale reference and anatomical orientation
- Files are saved in the "VoronoiWaveformMAPS" directory

## Configuration Options
- **Butter_Cut_Off**: Cutoff frequency for Butterworth filtering (default: 45Hz)
- **Order**: Filter order for Butterworth filter (default: 4)
- **Scale_lining**: Controls display of scale bar (1=on, 0=off)
- **Normalized_Train**: Toggles amplitude normalization (1=on, 0=off)
- **Stacked_filtration**: Enables post-averaging low-pass filtering (1=on, 0=off)
- **Pre_filt**: Enables pre-averaging high-pass filtering (1=on, 0=off)
- **Attinuation**: Controls scaling of inactive sites (1=attenuate, 0=normal scaling)
- **Scale_map_ratio**: Sets the scaling range for waveform amplitudes

## Signal Processing
- **Trimmed Mean**: Uses trimmed mean (40% trimming) for robust averaging
- **Filtering Options**:
  - Pre-filtering: Optional high-pass filtering at 20Hz
  - Post-filtering: Low-pass filtering at 45Hz
- **Amplitude Normalization**: Scales waveforms for visibility while maintaining relative magnitudes
- **Attenuation**: Optionally reduces amplitude of inactive sites for visual clarity

## Usage
1. Run the Voronoi_Waveform_plotter_V4.m script
2. Select the MetaData file when prompted
3. Select the SVG file for Voronoi cell definitions when prompted
4. The program will automatically process all channels for the selected amplitude
5. Output files will be saved in the "VoronoiWaveformMAPS" directory

## Technical Features
- Direct mapping of EMG waveforms to stimulation site coordinates
- Time alignment showing consistent stimulation periods
- Adaptive scaling based on Voronoi cell area calculations
- Two-pass drawing for layered visualization
- Sophisticated response validation logic for merged sites
- Activity markers distinguishing responsive from non-responsive sites

## Dependencies
- MATLAB (R2018b or later recommended)
- Signal Processing Toolbox (for filtering)
- Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform (for SVG processing)
- MetaData file from earlier processing stages

## Notes
- The visualization is most effective with filtered waveforms (Stacked_filtration=1)
- For publication-quality figures, SVG format is recommended
- Active sites are marked with "*-" before the site number
- Stimulation periods are highlighted with light gray rectangles
- The configuration options can be adjusted at the top of the script before running