# Cross-Case Heatmap (RDV22_CrossCASE_HeatmapV3.m)

## Overview
The Cross-Case Heatmap tool analyzes EMG patterns across multiple experimental subjects to identify conserved functional components in the sensorimotor cortex. By applying non-negative matrix factorization (NMF) to data from multiple cases, this tool reveals fundamental organizational principles that transcend individual anatomical variations.

## Requirements
- MATLAB R2019b or later
- Signal Processing Toolbox
- Statistics and Machine Learning Toolbox
- Wavelet Toolbox
- MiladFuncV1.m utility class
- Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform_CrossCase.m

## Input Data
The program requires the following data for each experimental case:
- Topo_Muscle_Cluster structure containing:
  - Signature_Refined: Waveform signatures for each cluster
  - Combined_Cluster_ident_Refined: Site identifiers for each cluster
  - Tracker_Refined: Cluster-to-channel mapping information
  - Combined_Verdict_Refined: Activity validation data
  - Case: Case identifier
  - Snips_fs: Sampling frequency
  - Size: Number of sites
  - Unblank_raw: Raw data configuration (0 or 1)

## Workflow

### 1. Data Collection and Integration
```matlab
% Batch upload of case data
for J=1:1:length(Case_set)
    % Load Topo files for current case
    [Partfile,Partpath] = uigetfile('*.mat',sprintf('load Topo files for %s',Case_set(J)));
    Load_Channel_name=sprintf('%s\\%s',Partpath,Partfile);
    load(Load_Channel_name);
    
    % Store case data in structured format
    fieldName = ['Case' num2str(J)];
    Superfile.(fieldName).Signatures = Topo_Muscle_Cluster.Signature_Refined;
    Superfile.(fieldName).Combined_Cluster_ident = Topo_Muscle_Cluster.Combined_Cluster_ident_Refined;
    % ... additional data storage
    
    % Build unified data matrices
    Super_Signatures = [Super_Signatures; Superfile.(fieldName).Signatures{i}];
    Super_Case_Tracker = [Super_Case_Tracker; (ones(Counter,1) * J)];
    % ... additional tracking variables
end
```

### 2. Signal Processing and Feature Extraction
```matlab
% Signal normalization and preprocessing
Testing_Channel_means = Super_Signatures-mean(Super_Signatures,2);
if Unblank_raw == 1 
    Testing_Channel_means = MiladFuncV1.Blank_it_Train(Testing_Channel_means);
end

% Segment signal into relevant windows
Testing_set_background = Testing_set(:,1:1526);
Testing_set_stim = Testing_set(:,1527:3053);
Testing_set_Post_stim = Testing_set(:,3054:end);
Testing_set_POST_ONSET = [Testing_set_stim, Testing_set_Post_stim];

% Prepare target signals for analysis
Target = normalize(Testing_set, 2, "range", [0 1]);
Target = Target.*Logic_verdict;
[activeIndices, activeSignals] = MiladFuncV1.separateActiveSignals(Target);
```

### 3. Non-negative Matrix Factorization (NMF)
```matlab
% Determine optimal number of NMF components
maxFactors = 20;
optimalFactors = MiladFuncV1.nnmf_Nfactor_error(NMF_Target, maxFactors);

% Perform NMF decomposition
[Low_order_approximation, SuperLeft, SuperRight] = MiladFuncV1.Milad_NMF(NMF_Target, Factors, algorithm);

% Calculate factor scores and case contributions
for i=1:1:size(SuperRight,1)
    Factor_Score(i,1) = sum(SuperLeft(:,i));
    for J=unique(Reduced_Super_Case_Tracker)'
        Case_Ratio(i,J) = sum(SuperLeft(Reduced_Super_Case_Tracker==J,i))/sum(SuperLeft(:,i));
    end
end
```

### 4. Case-Specific Analysis
```matlab
for Cas=unique(Reduced_Super_Case_Tracker)'
    % Extract case-specific data
    Locator = Reduced_Super_Case_Tracker == Cas;
    Left = SuperLeft(Locator>0,:);
    Right = SuperRight;
    
    % Calculate channel and site scores
    Norm_Left = normalize(Left, 1, 'range', [0 1]);
    for H=1:1:size(Norm_Left,2)
        Scores = Norm_Left(:,H);
        % Channel scores
        for i=1:1:16
            Channel_Score(i,H) = sum(Scores(Sub_Reduced_Tracker(:,1)==i));
        end
        % Site scores
        for i=1:1:Number_of_Sites
            % ... site score calculation
            Site_Score(i,H) = sum(Scores(mapper>0));
        end
    end
    
    % Generate Voronoi visualizations
    % ... visualization code
end
```

## Output
The program generates a comprehensive PDF report containing:
- Component waveform plots showing the characteristic shape of each NMF component
- Spatial heatmaps displaying component distribution across the cortical surface for each case
- Combined visualizations correlating spatial patterns with signal characteristics

## Parameters
Key parameters that may require adjustment:
- `Unblank_raw`: Set to 1 for raw data processing, 0 for preprocessed data
- `maxFactors`: Maximum number of components to consider in NMF analysis (default: 20)
- `algorithm`: NMF algorithm type ("mult" or "als")

## Usage Tips
1. Ensure that all case data follows the same structural format
2. Verify that all Voronoi map files are available for visualization
3. The quality of cross-case analysis depends on consistent recording and preprocessing across cases
4. For large datasets, consider increasing MATLAB's memory allocation
5. The optimal number of components may vary based on experimental design; review the error plot to confirm selection

## Troubleshooting
- If "Warning: Matrix is close to singular or badly scaled" appears, check normalization procedures
- If no patterns emerge across cases, verify that recording protocols were consistent
- For memory issues, consider processing cases in smaller batches

## References
- See manuscript for detailed methodology and interpretation guidelines
