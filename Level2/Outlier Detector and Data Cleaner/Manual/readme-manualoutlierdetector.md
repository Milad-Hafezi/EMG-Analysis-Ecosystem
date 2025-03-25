# Manual Outlier Detector (Meta_Data_plotter_flagger_V9_App.m)

## Description
The Manual Outlier Detector provides an interactive graphical interface for researcher-guided assessment of EMG signal quality. It enables users to visually inspect, compare, and select valid EMG trials while excluding outliers, leveraging human pattern recognition capabilities to identify complex artifacts and anomalies that might be challenging for algorithmic approaches.

## Inputs
- MetaData file (.mat) from the MetaData Generator component
- (Optional) Previously saved selection file for continuing interrupted assessment sessions

## Outputs
- Comprehensive output file (.mat) in the "CleanedData" directory containing:
  - `Mean_Cell`: Mean waveforms of selected valid trials for each channel, site, and amplitude
  - `Median_Cell`: Median waveforms of selected valid trials
  - `Selected_Train_number_cell`: Record of which specific trials were selected as valid
  - `Selected_Train_Matrix`: Binary matrix indicating selected (1) or excluded (NaN) trials
  - Processing metadata including type ("Manual_Clean") and version information

## Key Features
1. **Interactive Visual Interface**:
   - Stacked, normalized EMG traces for direct comparison
   - Checkboxes for individual trial selection
   - Clear visualization of stimulation periods
   - Adaptive layout handling different numbers of trials

2. **Selection Tools**:
   - Individual trial checkboxes
   - "Select All" option for efficient processing
   - Preservation of previous selections when editing
   - Trial identification and numbering for reference

3. **Navigation and Assessment**:
   - Dropdown menus for site, amplitude, and channel selection
   - "Plot Mean" button to preview selected trial averages
   - Comparison of selected vs. total dataset means
   - Both mean and median calculations for robustness

4. **Data Management**:
   - Arming/disarming toggle to prevent accidental saves
   - Option to start fresh or edit existing selections
   - Visual feedback for successful saving
   - Consistent output format compatible with automatic detection

## Usage
1. Run the Meta_Data_plotter_flagger_V9_App.m script
2. Select the MetaData file (.mat) when prompted
3. Choose to start fresh (0) or edit existing selections (any other key)
4. If editing, select the previous selection file
5. Use dropdown menus to navigate to desired channel, site, and amplitude
6. Visually assess the stacked EMG traces
7. Select valid trials using individual checkboxes (or "Select All")
8. Preview selections using the "Plot Mean" button
9. Arm the system (toggle button) and save selections using "Save changes"
10. Continue to next site/amplitude/channel using dropdown menus

## Dependencies
- MATLAB (R2018b or later recommended)
- MATLAB Signal Processing Toolbox (for Butterworth filtering)
- MATLAB App Designer components

## Notes
- The system requires at least one trial to be selected before saving
- The "Armed" status prevents accidental saving of selections
- Produces output files compatible with the Automatic Outlier Detector
- Adaptive layout works best with screen resolution 1024x768 or higher
- For large datasets, trials may be arranged in multiple columns
- Butterworth filtering (45Hz cutoff) is applied for visualization purposes