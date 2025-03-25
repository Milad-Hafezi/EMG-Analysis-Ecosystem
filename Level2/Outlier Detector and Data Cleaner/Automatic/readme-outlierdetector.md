# Automatic Outlier Detector (Automatic_Outlier_Detector_V3.m)

## Description
The Automatic Outlier Detector identifies and filters out anomalous EMG trials using advanced signal processing and machine learning techniques. This algorithm transforms EMG signals into time-frequency feature representations and employs unsupervised clustering to objectively separate genuine responses from artifacts and noise, ensuring only valid trials contribute to subsequent analyses.

## Inputs
- MetaData file (.mat) from the MetaData Generator component
- Optional configuration settings for processing parameters

## Outputs
- Comprehensive output file (.mat) in the "AutoCleanedData" directory containing:
  - `Mean_Cell`: Mean waveforms of valid trials for each channel, site, and amplitude
  - `Median_Cell`: Median waveforms for each condition
  - `Selected_Train_number_cell`: Record of which specific trials were selected as valid
  - `Selected_Train_Matrix`: Binary matrix indicating valid (1) or invalid (NaN) trials
  - Processing metadata and algorithm parameters
- Optional PDF reports visualizing clustering results (when PDF_logic=1)

## Key Algorithms
1. **Time-Frequency Analysis**:
   - Short-Time Fourier Transform (STFT) via the `STFT_Feature` function
   - Converts EMG signals into multidimensional feature representations
   - Captures both temporal evolution and frequency content

2. **Unsupervised Clustering**:
   - K-means clustering through the `KCluster` function
   - Separates trials into two clusters (genuine responses vs. artifacts/noise)
   - Multiple clustering replicates for robust separation

3. **Objective Cluster Evaluation**:
   - Davies-Bouldin Index calculation via the `DB_Indexing` function
   - Determines which cluster contains meaningful signals
   - Combines inter-cluster and intra-cluster distance metrics
   - Provides mathematical basis for distinguishing signal from noise

4. **Artifact Reduction**:
   - Optional blanking via the `Blank_it` function
   - Configurable blanking window parameters
   - Replaces artifact periods with pre-artifact values

## Usage
1. Run the Automatic_Outlier_Detector_V3.m script
2. Select the MetaData file (.mat) when prompted
3. The program will automatically:
   - Process each channel, amplitude, and site combination
   - Perform time-frequency feature extraction
   - Cluster trials and identify the valid cluster
   - Generate mean and median waveforms from valid trials
   - Create a record of which trials were selected
4. Results will be saved in the "AutoCleanedData" directory

## Dependencies
- MATLAB (R2018b or later recommended)
- MATLAB Signal Processing Toolbox (for spectrogram and filtering functions)
- MATLAB Statistics and Machine Learning Toolbox (for clustering)

## Notes
- The algorithm automatically determines the number of clusters (set to 2)
- Davies-Bouldin Index provides an objective measure of cluster quality
- The PDF report option (PDF_logic=1) creates visualizations for each condition
- Processing is performed hierarchically across channels, amplitudes, and sites
- The output serves as a cleaned dataset for all subsequent analyses
