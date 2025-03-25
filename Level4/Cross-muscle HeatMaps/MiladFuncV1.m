 classdef MiladFuncV1
    % MiladFuncV1 - A collection of utility functions
    
    methods (Static)           
                            
                            function Revised_Report=plotClusteringTree(data,activeIndices)
                            % PLOTCLUSTERINGTREE Plots a hierarchical tree from the given data.
                            %
                            %   data is an N x 1 struct array with fields:
                            %      - data.ClusterID : The cluster index (numeric).
                            %      - data.Level     : The depth level (numeric).
                            %      - data.Members   : A numeric array of members in that cluster.
                            %
                            % This function:
                            %   1) Finds parent–child relationships based on child's members being
                            %      a subset of the parent's members at the next level (Level+1).
                            %   2) Builds a digraph and uses a layered layout with a top-down direction.
                            %   3) Issues a warning (instead of an error) if no edges are found.
                                %%%%%%%%%%%%%%%%%%%%%%%
                                % Revising the member's labels
                                Revised_Report=data;

                                for JJ=1:1: size(data,1)
                                    All_memebrs = data(JJ).Members;
                                      Revised_Report(JJ).Members=activeIndices(All_memebrs);
                                end
                                
                                

                                % ---------------------------------------------------------------------
                                % 1) Build edges with default subset logic: childMembers ⊆ parentMembers
                                % ---------------------------------------------------------------------
                                edgeCell = {};   % to store {parentNode, childNode} pairs
                                allLevels = [data.Level];
                                minLvl    = min(allLevels);
                                maxLvl    = max(allLevels);
                              
                                for lvl = minLvl : maxLvl
                                    % indices of clusters at this level
                                    idxThisLevel = find([data.Level] == lvl);
                                    % indices of clusters at the next level
                                    idxNextLevel = find([data.Level] == (lvl + 1));
                            
                                    for i = idxThisLevel
                                        childCluster  = data(i);
                                        childMembers  = childCluster.Members;
                            
                                        for j = idxNextLevel
                                            parentCluster = data(j);
                                            parentMembers = parentCluster.Members;
                            
                                            % Subset check: childMembers ⊆ parentMembers
                                            if all(ismember(childMembers, parentMembers))
                                                parentNode = sprintf('L%d_C%d', parentCluster.Level, parentCluster.ClusterID);
                                                childNode  = sprintf('L%d_C%d', childCluster.Level,  childCluster.ClusterID);
                            
                                                % Store edge (parent -> child)
                                                edgeCell(end+1,:) = { parentNode , childNode }; %#ok<AGROW>
                                            end
                                        end
                                    end
                                end
                            
                                % ---------------------------------------------------------------------
                                % 2) Create a digraph and plot (no error on empty)
                                % ---------------------------------------------------------------------
                                if isempty(edgeCell)
                                    warning('No edges found. Plotting an empty or disconnected graph.');
                                    % Build a graph with the unique nodes but no edges
                                    allNodes = arrayfun(@(x) sprintf('L%d_C%d', x.Level, x.ClusterID), ...
                                                        data, 'UniformOutput', false);
                                    G = digraph();
                                    G = addnode(G, unique(allNodes));
                                else
                                    G = digraph(edgeCell(:,1), edgeCell(:,2));
                                end
                            
                                % ---------------------------------------------------------------------
                                % 3) Plot the graph
                                % ---------------------------------------------------------------------
                                figure('Name','Clustering Tree','Color','w');
                                plot(G, 'Layout','layered', 'Direction','down');
                                title('Clustering Tree (No Level Flipping, No Error on Empty)');
                                axis off
                            end

                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Hierarchical Clustering with Recompute and Progress Reporting
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            function [clusterAssignments, progressReport] = hierarchical_linkage_with_recompute(Testing_set_POST_ONSET, maxDepth, minClusterSize, maxWarpingWindow)
                                N = size(Testing_set_POST_ONSET,1);
                                clusterAssignments = zeros(N,1);    
                                initialIndices = (1:N)'; 
                                currentClusterID = 1;
                                progressReport = []; % Initialize progress report
                                [clusterAssignments, ~, progressReport] = MiladFuncV1.recursivelySplit(Testing_set_POST_ONSET, initialIndices, clusterAssignments, currentClusterID, maxDepth, minClusterSize, maxWarpingWindow, progressReport);
                            end
                            
                            function [clusterAssignments, nextClusterID, progressReport] = recursivelySplit(Testing_set_POST_ONSET, indices, clusterAssignments, currentClusterID, maxDepth, minClusterSize, maxWarpingWindow, progressReport)
                                n = length(indices);
                            
                                % Base condition: If no more depth or cluster too small, assign a final cluster ID
                                if maxDepth == 0 || n < minClusterSize
                                    clusterAssignments(indices) = currentClusterID;
                                    % Log the final cluster state
                                    progressReport = [progressReport; struct('ClusterID', currentClusterID, 'Level', maxDepth, 'Members', indices')];
                                    nextClusterID = currentClusterID + 1;
                                    return; 
                                end
                            
                                % Extract the subset of EMGs for this cluster
                                subsetEMGs = Testing_set_POST_ONSET(indices, :);
                                
                                % Compute a new DTW distance matrix for these EMGs
                                subDistMatrix = MiladFuncV1.computeDTWMatrix(subsetEMGs, maxWarpingWindow);
                            
                                % Convert to condensed distance vector for linkage
                                condDist = squareform(subDistMatrix, 'tovector');
                                
                                % Perform linkage on this newly computed DTW distance matrix
                                Z = linkage(condDist, 'average');
                                
                                % Cluster into two clusters
                                c = cluster(Z, 'maxclust', 2);
                                
                                cluster1_indices = indices(c == 1);
                                cluster2_indices = indices(c == 2);
                                
                                % If splitting did not produce meaningful two subclusters, assign one cluster
                                if isempty(cluster1_indices) || isempty(cluster2_indices)
                                    clusterAssignments(indices) = currentClusterID;
                                    % Log the cluster state
                                    progressReport = [progressReport; struct('ClusterID', currentClusterID, 'Level', maxDepth, 'Members', indices')];
                                    nextClusterID = currentClusterID + 1;
                                    return; 
                                end
                            
                                % Log progress for this split
                                progressReport = [progressReport; struct('ClusterID', currentClusterID, 'Level', maxDepth, 'Members', indices')];
                            
                                % Recursively split the first subcluster
                                [clusterAssignments, nextClusterID, progressReport] = MiladFuncV1.recursivelySplit(Testing_set_POST_ONSET, cluster1_indices, clusterAssignments, currentClusterID, maxDepth - 1, minClusterSize, maxWarpingWindow, progressReport);
                                
                                % Recursively split the second subcluster, starting IDs from nextClusterID
                                [clusterAssignments, nextClusterID, progressReport] = MiladFuncV1.recursivelySplit(Testing_set_POST_ONSET, cluster2_indices, clusterAssignments, nextClusterID, maxDepth - 1, minClusterSize, maxWarpingWindow, progressReport);
                                
                                return;
                            end
                            
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Reassigning Small Clusters with Reporting
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            function [updatedClusters, reassignmentReport] = reassignSmallClusters(clusterAssignments, Testing_set_POST_ONSET, MemberLimit,maxWarpingWindow)
                            % Find cluster sizes
                            uniqueClusters = unique(clusterAssignments);
                            uniqueClusters = uniqueClusters(uniqueClusters > 0);
                            
                            clusterSizes = arrayfun(@(c) sum(clusterAssignments == c), uniqueClusters);
                            
                            % Identify small clusters (fewer than MemberLimit members)
                            smallClusters = uniqueClusters(clusterSizes < MemberLimit);
                            reassignmentReport = struct('SmallCluster', {}, 'Member', {}, 'ReassignedTo', {}, 'Distances', {});
                            
                            if isempty(smallClusters)
                                % No small clusters to reassign
                                updatedClusters = clusterAssignments;
                                return;
                            end
                            
                            % Identify large clusters
                            largeClusters = uniqueClusters(clusterSizes >= MemberLimit);
                            
                            % If there are no large clusters, then we can't reassign to anything.
                            if isempty(largeClusters)
                                updatedClusters = clusterAssignments;
                                return;
                            end
                            
                            % Compute the mean waveform of each large cluster
                            clusterMeans = cell(numel(largeClusters),1);
                            for i = 1:numel(largeClusters)
                                clID = largeClusters(i);
                                members = Testing_set_POST_ONSET(clusterAssignments == clID,:);
                                clusterMeans{i} = mean(members,1); % mean waveform
                            end
                            
                            % Now dissolve each small cluster
                            for sc = smallClusters'
                                % Get members of small cluster
                                scMembers = find(clusterAssignments == sc);
                                
                                for m = scMembers'
                                    % EMG waveform of this member
                                    emgSignal = Testing_set_POST_ONSET(m,:);
                                    
                                    % Compute DTW distance of this EMG to each large cluster mean
                                    dists = zeros(numel(largeClusters),1);
                                    for i = 1:numel(largeClusters)
                                        meanSignal = clusterMeans{i};
                                        dists(i) = dtw(emgSignal, meanSignal,maxWarpingWindow);

                                        
                                    end
                                    
                                    % Find closest large cluster
                                    [~, minIdx] = min(dists);
                                    closestCluster = largeClusters(minIdx);
                                    
                                    % Assign this EMG to the closest large cluster
                                    clusterAssignments(m) = closestCluster;
                                    
                                    % Log the reassignment
                                    reassignmentReport = [reassignmentReport; struct('SmallCluster', sc, 'Member', m, 'ReassignedTo', closestCluster, 'Distances', dists')];
                                end
                            end
                            
                            % After reassigning, some clusters (the small ones) are now empty.
                            % Remove empty clusters and renumber them from 1 to K
                            finalClusters = clusterAssignments;
                            uniqueFinal = unique(finalClusters);
                            uniqueFinal = uniqueFinal(uniqueFinal > 0);
                            
                            % Map old cluster IDs to new ones
                            newIDMap = containers.Map(num2cell(uniqueFinal), num2cell(1:numel(uniqueFinal)));
                            
                            % Ensure that only valid keys are used for remapping
                            updatedClusters = zeros(size(finalClusters));
                            for i = 1:length(finalClusters)
                                if isKey(newIDMap, finalClusters(i))
                                    updatedClusters(i) = newIDMap(finalClusters(i));
                                else
                                    % Handle unexpected cases where a cluster ID is missing in newIDMap
                                    updatedClusters(i) = 0; % Assign to a default cluster or handle appropriately
                                end
                            end
                        end

        
                                
                                
                     
        
        
        
        
                                function Square_dist_dtw = computeDTWMatrix(focused_signal, maxWarpingWindow)
                                % computeDTWMatrixParallel - Computes the DTW distance matrix for EMG signals using parallel processing.
                                %
                                % Inputs:
                                %   focused_signal - A matrix where each row is an EMG signal.
                                %   maxWarpingWindow - Maximum warping window for DTW (e.g., 150 samples).
                                %
                                % Output:
                                %   Square_dist_dtw - Symmetric matrix with DTW distances between signals.
                            
                                % Number of signals
                                N = size(focused_signal, 1);
                            
                                % Precompute pair indices
                                [rowIdx, colIdx] = find(triu(true(N), 1)); % Upper triangle indices
                                numPairs = length(rowIdx);
                            
                                % Preallocate distance vector
                                dist_dtw = zeros(numPairs, 1);
                            
                                % Parallel computation
                                parfor k = 1:numPairs
                                    i = rowIdx(k);
                                    j = colIdx(k);
                                    % Compute DTW distance for pair (i, j)
                                    dist_dtw(k) = dtw(focused_signal(i, :), focused_signal(j, :), maxWarpingWindow);
                                end
                            
                                % Convert to symmetric matrix
                                Square_dist_dtw = zeros(N, N);
                                for k = 1:numPairs
                                    i = rowIdx(k);
                                    j = colIdx(k);
                                    % Populate symmetric matrix
                                    Square_dist_dtw(i, j) = dist_dtw(k);
                                    Square_dist_dtw(j, i) = dist_dtw(k);
                                end
                            end


                        function [activeIndices, activeSignals] = separateActiveSignals(Target)
                            % Separate active and negative EMG signals
                            Logic_verdict = (mean(Target, 2) > 0); % Logic to separate active (1) and negative (0)
                            activeIndices = find(Logic_verdict); % Indices of active EMG signals
                            activeSignals = Target(activeIndices, :); % Active EMG signals themselves
                        end
        
        
        


                    
                    function featureTable = extractEMGFeatures(emgMatrix, samplingFrequency)
                    % extractEMGFeatures - Extracts time-domain and frequency-domain features
                    % from an EMG signal matrix where each row is an EMG signal.
                    %
                    % Inputs:
                    %   emgMatrix - A matrix where each row is an EMG signal.
                    %   samplingFrequency - The sampling frequency of the EMG signals (Hz).
                    %
                    % Output:
                    %   featureTable - A table with rows representing signals and columns as features.
                    
                    % Number of signals
                    numSignals = size(emgMatrix, 1);
                    
                    % Initialize feature matrix
                    featureMatrix = zeros(numSignals, 9); % 5 time-domain + 4 frequency-domain features
                    
                    % Feature extraction
                    for i = 1:numSignals
                        signal = emgMatrix(i, :);
                        
                        % Time-domain features
                        lyapunovExponent = MiladFuncV1.calculateLyapunovExponent(signal);
                        correlationDimension = MiladFuncV1.calculateCorrelationDimension(signal);
                        approximateEntropy = MiladFuncV1.calculateApproxEntropy(signal);
                        shapeFactor = rms(signal) / mean(abs(signal));
                        impulseFactor = max(abs(signal)) / mean(abs(signal));
                        
                        % Frequency-domain features (using autoregressive model)
                        [frequencies, bandPower] = MiladFuncV1.calculateFrequencyFeatures(signal, samplingFrequency);
                        
                        % Store features
                        featureMatrix(i, :) = [lyapunovExponent, correlationDimension, approximateEntropy, ...
                                               shapeFactor, impulseFactor, frequencies, bandPower];
                    end
                    
                    % Create labeled table
                    featureTable = array2table(featureMatrix, ...
                        'VariableNames', {'LyapunovExponent', 'CorrelationDimension', ...
                                          'ApproxEntropy', 'ShapeFactor', 'ImpulseFactor', ...
                                          'Freq1', 'Freq2', 'Freq3', 'BandPower'});
                end
                
                % Subfunction: Calculate Lyapunov Exponent
                function lyapExp = calculateLyapunovExponent(signal)
                    % Placeholder for Lyapunov Exponent calculation
                    % Implement proper Lyapunov Exponent calculation based on your requirements
                    lyapExp = rand(); % Replace with actual computation
                end
                
                % Subfunction: Calculate Correlation Dimension
                function corrDim = calculateCorrelationDimension(signal)
                    % Placeholder for Correlation Dimension calculation
                    % Implement proper Correlation Dimension calculation based on your requirements
                    corrDim = rand(); % Replace with actual computation
                end
                
                % Subfunction: Calculate Approximate Entropy
                function approxEntropy = calculateApproxEntropy(signal)
                    % Placeholder for Approximate Entropy calculation
                    % Implement proper Approximate Entropy calculation based on your requirements
                    approxEntropy = rand(); % Replace with actual computation
                end
                
                % Subfunction: Calculate Frequency-Domain Features
                function [frequencies, bandPower] = calculateFrequencyFeatures(signal, fs)
                    % Estimate the power spectral density using an autoregressive model
                    order = 8; % AR model order
                    [arCoeffs, noiseVariance] = aryule(signal, order);
                    [psd, freqs] = freqz(sqrt(noiseVariance), arCoeffs, 1024, fs);
                    
                    % Find the first three natural frequencies (peaks in the PSD)
                    [~, peakIdx] = findpeaks(abs(psd), 'SortStr', 'descend');
                    naturalFrequencies = freqs(peakIdx(1:min(3, length(peakIdx)))); % Take up to 3 peaks
                    
                    % Calculate band power (total power in the spectrum)
                    bandPower = bandpower(real(psd), freqs, 'psd');
                    
                    % Pad natural frequencies if fewer than 3 are found
                    frequencies = zeros(1, 3);
                    frequencies(1:length(naturalFrequencies)) = naturalFrequencies;
                end



           

                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            % Cluster based on Lags
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                            function updatedCluster = reassignByAvgLag(activeSignals, Reduced_cluster,Latency_vector_active, maxLag)
%REASSIGNBYAVGLAG Reassign signals to new clusters if their average lag 
% (w.r.t. other signals in the same cluster) is too large.
%
%   updatedCluster = reassignByAvgLag(activeSignals, reducedCluster, maxLag)
%
%   INPUTS:
%       activeSignals  - N x M matrix of signals (N=number of signals)
%       reducedCluster - N x 1 vector of initial cluster assignments
%       maxLag         - scalar threshold for allowable average signal delay
%
%   OUTPUT:
%       updatedCluster - N x 1 vector of new cluster labels after 
%                       reassigning signals whose average lag is either
%                       above maxLag or significantly larger than the 
%                       cluster's mean average-lag.
%
%   NOTE:
%       This function only processes the *original* clusters in reducedCluster;
%       newly created clusters are not further split in the same call.

    % Copy the input cluster labels for output
    updatedCluster = Reduced_cluster;
    
    % Identify only the original clusters
    originalClusters = unique(Reduced_cluster, 'stable');
    
    % Define a factor for deciding “significantly different”
    % e.g., if average lag > factor * (cluster avg-lag) => outlier
%    factor = 1.5;   % <-- adjust to your needs
    
    for cIdx = 1:numel(originalClusters)
        currentClusterID = originalClusters(cIdx);
        
        % Find members of this cluster (in the updated assignment so far)
        memberIdx = find(updatedCluster == currentClusterID);
        
        if numel(memberIdx) < 2
            % If only one member, no pairwise comparison needed
            continue;
        end
        
        % ------------------------------------------------------------
        % 1) Compute average-lag for each signal in this cluster
        % ------------------------------------------------------------
        
        % Preallocate
        avgLagAll = zeros(size(memberIdx));
        clear Abs_latency MaxLagAll
        for m = 1:numel(memberIdx)
            sigM = memberIdx(m);
            
            % Compute the lag of sigM vs all others in the cluster
            lagSum = 0;
            count  = 0;
            clear delays 
            for n = 1:numel(memberIdx)
                if n == m; continue; end
                sigN = memberIdx(n);
                
                delayMN =  MiladFuncV1.computeSignalDelay(activeSignals(sigM, :), ...
                                             activeSignals(sigN, :));
                lagSum = lagSum + delayMN;
                count  = count + 1;
                delays(n)=delayMN;

            end
            
            % Average lag for this signal
            avgLagAll(m) = lagSum / count;
            MaxLagAll(m)=max(delays);
            Abs_latency(m)=Latency_vector_active(sigM);
            
        end
        
        % ------------------------------------------------------------
        % 2) Compute the cluster’s mean average-lag
        % ------------------------------------------------------------
        clusterMeanAvgLag = mean(avgLagAll);
%            figure
%          for ii = 1:numel(memberIdx)
%               sigN = memberIdx(ii);
%              plot(activeSignals(sigN, :), 'DisplayName',num2str(sigN));
%              hold on
%          end
%          legend

        % ------------------------------------------------------------
        % 3) Check each signal for outlier conditions
        % ------------------------------------------------------------
        for m = 1:numel(memberIdx)
            sigM = memberIdx(m);
            
            % If already reassigned in some earlier iteration, skip it
            if updatedCluster(sigM) ~= currentClusterID
                continue;
            end
            
            sigAvgLag = avgLagAll(m);
            max_lag=MaxLagAll(m);
            
            % Condition A: average lag is above maxLag
            condA = (max_lag > maxLag);
            % Condition B: average lag is “significantly different” from
            %             the cluster mean (e.g., 1.5 times bigger)
                
%                 condB = (sigAvgLag > factor * clusterMeanAvgLag);  % I disabled condB because it was not robustly perforing
                     

            % Condition C: outlier detection
                 Outliers_Avg=isoutlier(avgLagAll);
                 condC=Outliers_Avg(m);
                 % condition D:
             if numel(memberIdx)==2
                  condD= diff(Abs_latency) > (maxLag/4);
                  ForceMajor=1;

             else
                  Outliers_Lat=isoutlier(Abs_latency);
                  condD=Outliers_Lat(m);
                  ForceMajor=0;
             end

%               (condA) && (condC || condD)
% 
%             (condD) && (condC || condA)

%            (condD) || (condC || condA)
%            
%            ((condA) && (condC || condD))   || (ForceMajor && condD)

            if ((condA) && (condC || condD))   || (ForceMajor && condD)
                % Create a new cluster ID
                newClusterID = max(updatedCluster) + 1;
                
                % Reassign this signal to the new cluster
                updatedCluster(sigM) = newClusterID;
                
            end
        end
    end
end

%% Example placeholder for the delay computation function
function delayVal = computeSignalDelay(signalA, signalB)
%COMPUTESIGNALDELAY Returns a scalar measure of delay between two signals.
% Replace this logic with your actual cross-correlation, time-shift, etc.

    % For demonstration, we measure delay as the lag at max cross-correlation:
    [r, lags] = xcorr(signalA, signalB,750);
    [~, maxIdx] = max(abs(r));
    delayVal = abs(lags(maxIdx));
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_latency(Testing_set,Latencies)
                non_nan=~isnan(Latencies);

                
            for i=1:1:length(Latencies)
                if non_nan(i)==1
                    figure
                    plot(Testing_set(i,:),'DisplayName',num2str(i))
                    hold on
                    xline(Latencies(i),'r',num2str(Latencies(i)))
                    legend
                end
            end
        end
              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Onset-Detection Utility: Baseline Statistics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [meanVals, stdVals] = getBaselineStats(Testing_set, baselineCols)
            % GETBASELINESTATS Compute mean and std for each row in the given baseline region.
            %   [meanVals, stdVals] = MiladFuncV1.getBaselineStats(Testing_set, baselineCols)
            %
            %   Testing_set: N x T matrix (N signals, T samples)
            %   baselineCols: vector of column indices for baseline region (e.g. 1:1526)

            baselineData = Testing_set(:, baselineCols);
            meanVals = mean(baselineData, 2);  % N x 1
            stdVals  = std(baselineData, 0, 2);% N x 1
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1) Double-Threshold Method
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function onsetSamples = doubleThresholdOnsetAll(Testing_set, lowFactor, highFactor, baselineCols, minPersistSamples)
            % DOUBLETHRESHOLDONSETALL Detect onsets using a double-threshold approach.
            %
            %   onsetSamples = MiladFuncV1.doubleThresholdOnsetAll(Testing_set, ...
            %                      lowFactor, highFactor, baselineCols, minPersistSamples)
            %
            %   - lowFactor, highFactor: multipliers for baseline STD
            %   - baselineCols: columns where signal is quiet (e.g. 1:1526)
            %   - minPersistSamples: # consecutive samples above the low threshold 
            %                        required before confirming onset
            %   - onsetSamples: 1 x N array of onset indices (NaN if none found)

            [meanVals, stdVals] = MiladFuncV1.getBaselineStats(Testing_set, baselineCols);
            [numSignals, numSamples] = size(Testing_set);

            onsetSamples = nan(1, numSignals);
            
            for iSig = 1:numSignals
                sig = Testing_set(iSig, :);
                
                lowThresh  = meanVals(iSig) + lowFactor  * stdVals(iSig);
                highThresh = meanVals(iSig) + highFactor * stdVals(iSig);
                foundOnset = false;
                
                for n = 1:numSamples - minPersistSamples
                    if all(sig(n:n+minPersistSamples-1) > lowThresh)
                        % Potential onset region
                        % Check if within next ~50 samples we exceed the high threshold
                        searchWindowEnd = min(n+50, numSamples);
                        if any(sig(n:searchWindowEnd) > highThresh)
                            onsetSamples(iSig) = n;
                            foundOnset = true;
                            break;
                        end
                    end
                end
                
                if ~foundOnset
                    onsetSamples(iSig) = NaN;
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4) Adaptive Threshold (Mean + k*Std) with Min Duration
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function onsetSamples = adaptiveThresholdAll(Testing_set, baselineCols, factor, minDuration)
            % ADAPTIVETHRESHOLDALL Detect onsets using baseline mean + factor*std.
            %
            %   onsetSamples = MiladFuncV1.adaptiveThresholdAll(Testing_set, ...
            %                       baselineCols, factor, minDuration)

            [meanVals, stdVals] = MiladFuncV1.getBaselineStats(Testing_set, baselineCols);
            [numSignals, numSamples] = size(Testing_set);

            onsetSamples = nan(1, numSignals);

            for iSig = 1:numSignals
                sig = Testing_set(iSig, :);
                
                thresh = meanVals(iSig) + factor * stdVals(iSig);

                for n = 1:numSamples - minDuration
                    if all(sig(n:n+minDuration-1) > thresh)
                        onsetSamples(iSig) = n;
                        break;
                    end
                end
            end
        end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 6) Wavelet-Based Approach (Example)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function onsetSamples = waveletOnsetAll(Testing_set, baselineCols, waveletName, level, factor, minDuration)
            % WAVELETONSETALL Detect onsets by analyzing wavelet detail coefficients.
            %
            %   onsetSamples = MiladFuncV1.waveletOnsetAll(Testing_set, ...
            %                      baselineCols, waveletName, level, factor, minDuration)

            [numSignals, numSamples] = size(Testing_set);
            onsetSamples = nan(1, numSignals);

            for iSig = 1:numSignals
                sig = Testing_set(iSig, :);

                % Wavelet decomposition
                [C, L] = wavedec(sig, level, waveletName);
                
                % Reconstruct detail signals from levels 2 up to 'level'
                bandSig = zeros(1, numSamples);
                for lvl = 2:level
                    D = wrcoef('d', C, L, waveletName, lvl);
                    bandSig = bandSig + D;
                end

                % Rectify if needed
                bandSig = abs(bandSig);

                % Baseline stats
                baselineVals = bandSig(baselineCols);
                meanB = mean(baselineVals);
                stdB  = std(baselineVals);

                threshold = meanB + factor * stdB;

                % Find first crossing that persists
                for n = 1:numSamples - minDuration
                    if all(bandSig(n:n+minDuration-1) > threshold)
                        onsetSamples(iSig) = n;
                        break;
                    end
                end
            end
        end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            function [Signatures,Tracker,Combined_Cluster_ident,Combined_Verdict]=properCluster(MyTarget,Cluster_Info,Verdictor)
                            
                            
                                            clear Signatures Cluster_ident
                                            for J=1:1:size(MyTarget,2) %channels
                                                clear Trends Memmer_IDs
                                                for k=1:1:size(MyTarget(:,1),1)
                                                    Members=MyTarget{k,J};
                                                    if ~isempty(Members)
                                                        Trends(k,:)= trimmean(Members,30,1);
                                                        Memmer_IDs{k,:}=Cluster_Info{k,J};
                                                    end
                                                end
                                              
                                                Signatures{1,J}=Trends;
                                                Cluster_ident{1,J}=Memmer_IDs;
                                            end
                                            
                                            % by this point, we should have a cell with all trends of a signal muscle
                                            % building the combined Matrix
                                            
                                            % Initialize combined matrix and a tracker matrix
                                            Combined_Muscles = [];
                                            Tracker          = [];
                                            for i = 1:length(Signatures)
                                                currentBlock = Signatures{i}; 
                                                numRows      = size(currentBlock,1);
                                                
                                                Combined_Muscles = [Combined_Muscles; currentBlock];
                                                
                                                % Column 1 = cell index (i)
                                                % Column 2 = row index within the i-th cell
                                                rowIndices = (1:numRows).';  % row indices in current cell
                                                Tracker    = [Tracker; [repmat(i,numRows,1), rowIndices]];
                                            
                                            end
                                            
                                            clear Combined_Cluster_ident
                                            for no=1:1:size(Tracker,1)
                                                   clear Sites
                                                   Sites=Cluster_ident{Tracker(no,1)}{Tracker(no,2), 1};
                                                   Combined_Cluster_ident{no,1}=Sites;
                                                   Combined_Verdict(no,1)=sum(Verdictor(Sites,Tracker(no,1)));
                                            end
                            end




            
            function Struct=Get_Train(Selected_Channel,god,Amp_interest)
            global Meta_Data Unblank_raw
            
            Guide_trains  = Meta_Data{1, Selected_Channel}.Guide_trains ; 
            Guide_Matrix  = Meta_Data{1, Selected_Channel}.Guide_Matrix;
            % Muscle_name   = Meta_Data{1, Selected_Channel}.Muscle_name;
            %       Case    = Meta_Data{1, Selected_Channel}.Case;
            % Channel_Number= Meta_Data{1, Selected_Channel}.Channel_Number;
               Snips_fs    = Meta_Data{1, Selected_Channel}.Snips_fs;
            
                                 [row,col, val] =find (Guide_Matrix(:,6)==Amp_interest);
                                                  
            %                      clear Sub_matrix_Trains Sub_matrix_info
                                 Sub_matrix_info= Guide_Matrix(row,:);
                                 Sub_matrix_Trains=Guide_trains(row,:);
            
                                        %% Amplitude loop
                                        
                                                % finding the correct sites
                                   [ro,co,v]=find (Sub_matrix_info(:,1)==god);
                                    
                                            if isempty(ro)
                                                disp('Amp not found for this site')
                                            else
                                                Train_info=Sub_matrix_info(ro,:);
                                                Trains=Sub_matrix_Trains(ro,:);
                                            end
                          Trace=sprintf('Ch%g-Site%g-Amp%g-Raw%g',Selected_Channel,god,Amp_interest,Unblank_raw);
                          Struct.Train_info=Train_info;
                          Struct.Trains=Trains;
                          Struct.Trace=Trace;
                          Struct.Precapture_index=Meta_Data{1, Selected_Channel}.Guide_Matrix(1,9);
                          Struct.Snips_fs= Snips_fs;
            
            
            
            end
            
            
            
            
            
                  
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            
            
            function Stack_plotter(Struct)
                % Check if the input is a valid matrix
                matrix=Struct.Trains;
                if ~ismatrix(matrix)
                    error('Input must be a matrix.');
                end
                
                % Get the size of the matrix
                [m, n] = size(matrix);
                
                % Create a new figure
                figure;
                
                % Loop through each row
                for i = 1:m
                    % Get the current row
                    row = matrix(i, :);
                    
                    % Calculate the offset based on the maximum of the current row
                    offset = max(matrix(:)) * (m - i);
                    
                    % Plot the current row, adding the offset
                    plot(row + offset, 'LineWidth', 1);
                    hold on; % Hold on to plot multiple lines
                    
                    % Optional: Add labels and title
                    xlabel('Column Index');
                    ylabel('Value');
                    title(sprintf('%s',Struct.Trace));
                    
                    % Add a grid for better visualization
                    grid on;
                end
                
                % Adjust y-axis limits to fit all rows nicely
                ylim([0, max(matrix(:)) * m]);
                
                % Create a legend for clarity
                %legend(arrayfun(@(x) sprintf('Row %d', x), 1:m, 'UniformOutput', false));
                legend off
                grid off
                axis off
                xline(1526)
                xline(3052)
                hold off; % Release the hold
            end
            
            
            
            
            
            
            function  [U,S,V]=PCA_SVD(Trains,Trace)
            
                    obs=Trains;
                    titles=Trace;
            
                    [U,S,V]=svd(obs, 'econ');
                    
            %         figure
            %         
            %         
            %         subplot(1,2,1)
            %         semilogy(diag(S),'k-o' , 'LineWidth' , 2.5)
            %         set(gca, 'FontSize' , 15), axis tight, grid on
            %         subplot(1,2,2)
            %         
            %         plot (cumsum(diag (S)) , ' k-o' , 'LineWidth' , 2.5)
            %         set(gca, 'FontSize' , 15), axis tight, grid on
            %         set(gcf, 'Position', [1400 100 3*600 3*250] )
                    
                    figure,
                    
                    subplot(3,2,[2,4,6])
                    for i=1:size(obs,1)
                    
                            x=V(:,1)'*obs(i,:)';
                            y=V(:,2)'*obs(i,:)';
                            z=V(:,3)'*obs(i,:)';
                    
                            plot3(x,y,z,'rx','LineWidth',3)
                            view([22.670 -13.162])
                            hold on
            
                            X(i)=x;
                            Y(i)=y;
                            Z(i)=z;
                            
            
                    end
                    
                    xlabel('PC1'),ylabel('PC2'), zlabel('PC3')
                    title(sprintf('%s_ %g data points',titles,size(obs,1)))
                    grid on
            
            
                    
                      subplot(3,2,1)
                    plot(X)
                    hold on
            %         plot(abs(X))
                    title('PC1')
                     xline(1526)
                     xline(3052)
                     subplot(3,2,3)
                    plot(Y)
                     hold on
            %         plot(abs(Y))
                    title('PC2')
                    xline(1526)
                    xline(3052)
                     subplot(3,2,5)
                    plot(Z)
                    hold on
            %         plot(abs(Z))
                    title('PC3')
                    xline(1526)
                    xline(3052)
            
            end
            
            
            
            
            function [Swap_Struct,Blanked_Train,Blank_Throughput] =Blank_it(Struct)
                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
                                         Pulse_Period=4.99712;   %ms     0.4897  0.00499712
                                         Blanking_window=3; % ms
                                         Blanking_Shift=0; % a number from 0 to 100 (percentage) this can shift all the blanking windows by a fraction of the pulse period length. 
                                         forwarding=0;
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
            
                        Trains=Struct.Trains;
                        Precapture_index=Struct.Precapture_index;         
                        Snips_fs=Struct.Snips_fs;
            
            
                         % Logic vector:
                        
                        % Blank_Logic= zeros(1,size(Train,2));
                        % finding zero index
                        Blanking_onset_index=Precapture_index; %ceil
                        Blank_Steps=(Blanking_window/1000)*Snips_fs; %floor
                        Blank_Throughput=ceil(((Pulse_Period-Blanking_window)/1000)*Snips_fs);
                        Window=(Pulse_Period/1000)*Snips_fs;
                        
                        Blanking_Shift=(Blanking_Shift/100)*Window; %redefining the blanking_shift in index domain
                        Blanked_Train=Trains;
                        
                        for bil=1:1:100 % make sure to change this using Count matrix later if you are doing a short train stim
                                Starting=floor(Blanking_Shift +Blanking_onset_index+forwarding); 
                                Ending= ceil(Blanking_Shift+Blanking_onset_index+forwarding+Blank_Steps);
                                    for rows=1:1:size(Trains,1)
                                      Blanked_Train(rows,Starting:Ending)=Blanked_Train(rows,Starting-1);
                                    end
                        
                        %         Blank_Logic (Starting:Ending)=1;
                                forwarding=round(bil*Window);
                        end
                        % for pip=1:1:5
                        % figure
                        %  plot(Train(pip,:),'r')
                        %  hold on 
                        % plot(Blanked_Train(pip,:),'b')
                        % end
                        Swap_Struct=Struct;
                        Swap_Struct.Trains=Blanked_Train;
            
            
            
            end
            
            
            function Blanked_Train =Blank_it_Train(Trains)
            global Meta_Data Unblank_raw
                    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
                                         Pulse_Period=4.99712;   %ms     0.4897  0.00499712
                                         Blanking_window=3; % ms
                                         Blanking_Shift=0; % a number from 0 to 100 (percentage) this can shift all the blanking windows by a fraction of the pulse period length. 
                                         forwarding=0;
                          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                        
            
            %             Trains=Struct.Trains;
                        Precapture_index=Meta_Data{1, 1}.Guide_Matrix(1,9) ;
                        Snips_fs    = Meta_Data{1, 1}.Snips_fs; 
                      
            
            
                         % Logic vector:
                        
                        % Blank_Logic= zeros(1,size(Train,2));
                        % finding zero index
                        Blanking_onset_index=Precapture_index; %ceil
                        Blank_Steps=(Blanking_window/1000)*Snips_fs; %floor
                        Blank_Throughput=ceil(((Pulse_Period-Blanking_window)/1000)*Snips_fs);
                        Window=(Pulse_Period/1000)*Snips_fs;
                        
                        Blanking_Shift=(Blanking_Shift/100)*Window; %redefining the blanking_shift in index domain
                        Blanked_Train=Trains;
                        
                        for bil=1:1:100 % make sure to change this using Count matrix later if you are doing a short train stim
                                Starting=floor(Blanking_Shift +Blanking_onset_index+forwarding); 
                                Ending= ceil(Blanking_Shift+Blanking_onset_index+forwarding+Blank_Steps);
                                    for rows=1:1:size(Trains,1)
                                      Blanked_Train(rows,Starting:Ending)=Blanked_Train(rows,Starting-1);
                                    end
                        
                        %         Blank_Logic (Starting:Ending)=1;
                                forwarding=round(bil*Window);
                        end
                        % for pip=1:1:5
                        % figure
                        %  plot(Train(pip,:),'r')
                        %  hold on 
                        % plot(Blanked_Train(pip,:),'b')
                        % end
                       
                
            end
            
            
            function clusters=KClust(Y)
            
            
            
            if size(Y,1)>9
                ender=9;
            else
                ender=size(Y,1);
            end
            
            clust=zeros(size(Y,1),ender);
            
            for K=1:ender
            clust(:,K)= kmeans(Y, K);
            end
            eva=evalclusters(Y,clust,"CalinskiHarabasz");
            figure
            plot(eva)
            title('Optimat number of clusters K');
            
            clusters=kmeans(Y,eva.OptimalK,"Start","sample");
            figure
            
            scatter3(Y(:,1),Y(:,2),Y(:,3),50,clusters,'filled')
            title(sprintf('Kmeans clusters n=%g',eva.OptimalK));
            view(-50,8)
            
            end
            
            
            function  idx=BDSCANclust(Y)
            idx=dbscan(Y,10,10);
            
            figure
            scatter3(Y(:,1),Y(:,2),Y(:,3),50,idx,'filled')
            view(-50,8)
            title(sprintf('BDSCAN clusters n=%g',size(unique(idx),1)));
            
            end
            
            
            function Filtered_train=filtera(Trains,Low_pass_cutoff)
            global Meta_Data Unblank_raw
            Snips_fs    = Meta_Data{1, 1}.Snips_fs;
            
            Order=4;
                % Design the lowpass filter
                Wn=Low_pass_cutoff/(Snips_fs/2); % Filter parameters
                        [b,a] = butter(Order,Wn,'low');
            
                 for j=1:1:size(Trains,1)
            
                        Butter_Train(j,:)= filtfilt(b,a,Trains(j,:));
            
                  end
                    if exist('Butter_Train')
                            Filtered_train= Butter_Train;
                    end
            
            end
            
            
            
            function optimalFactors=nnmf_Nfactor_error(Testing_Channel_means,maxFactors)
                    
                    %% finding the optimal NMF factor
            %         The code determines the optimal number of factors for NNMF by splitting the data into training and validation sets, then averaging reconstruction errors over multiple runs for each factor count. It plots these errors and identifies the number of factors that yields the lowest reconstruction error, ensuring robust model evaluation.
                    
                    % Define the range of factors to test
            %         maxFactors; % Adjust based on your data
                    errors = zeros(maxFactors, 1);
                    
                    % Loop over different numbers of factors
                    for Factors = 1:maxFactors
                        [Left, Right] = nnmf(Testing_Channel_means, Factors);
                        
                        % Calculate reconstruction
                        Low_order_approximation = Left * Right;
                        
                        % Calculate reconstruction error (Frobenius norm)
                        errors(Factors) = norm(Testing_Channel_means - Low_order_approximation, 'fro');
                    end
                    
                    % % Plot the reconstruction error vs. number of factors
                    % figure;
                    % plot(1:maxFactors, errors, '-o');
                    % xlabel('Number of Factors');
                    % ylabel('Reconstruction Error');
                    % title('Reconstruction Error vs. Number of Factors');
                    % grid on;
                    
                    % Optionally, find the optimal number of factors
                    [~, optimalFactors] = min(errors); % Finds the index (number of factors) with the minimum error
                    
            end
            
            function [Low_order_approximation,Left,Right]=Milad_NMF(Testing_Channel_means,Factors,algorithm)
            
                    rng(1) % this will cause the reprocuability. Otherwise NMF would generate random combinations of signals.
                    [Left,Right]=nnmf(Testing_Channel_means,Factors,"algorithm",algorithm);
                    
                    figure
                    for j=1:1:size(Right,1)
                            if j==1
                            a=j;
                            b=j+1;
                            else
                            a=a+2;
                            b=b+2;
                            end
                            subplot(size(Right,1),2,a)
                            plot(Right(j,:))
                        
                            subplot(size(Right,1),2,b)
                            plot(normalize(Left(:,j),1,"range"))
                    end
                    
                    Low_order_approximation=Left*Right; % looks like the reduced order, has less noise and could be reused
            
            
            end
            
            
            
            
            
            
            function p=plot_voronoi()
            global Extract
            
            %% Color pallet
                Grey_Back=  [0.909803921568627   0.905882352941176   0.898039215686275];
                Grey_lines=  [0.819607843137255   0.827450980392157   0.831372549019608];
                Grey_Unfair= [0.498039215686275   0.498039215686275   0.498039215686275];
                
                        % plot the voronois
                             fig = figure;
                                     hold on;
                                     for id = 1:size(Extract,1)
                                        plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                                        p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                                        hold on;
                                        p.EdgeColor=Grey_lines;
                                       %plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                                       
                                       %text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
                                                            
                                     end
                                  
                                    axis off
                                    axis equal
                                        
            %                         xl = xlim; % getting the X and Y limits to put in for later plots
            %                         yl= ylim;
                        
                                    hold on
            
            
            end
            
            function p=Subplot_voronoi(Site_Channel_map,Cluster_members,fontsize)
            global Extract
            
            %% Color pallet
                Grey_Back=  [0.909803921568627   0.905882352941176   0.898039215686275];
                Grey_lines=  [0.819607843137255   0.827450980392157   0.831372549019608];
                Grey_Unfair= [0.498039215686275   0.498039215686275   0.498039215686275];
                
                        % plot the voronois
%                              fig = figure;
                                     hold on;
%                                      subplot(2,1,1)
                                      hold on;
                                     for id = 1:size(Extract,1)
                                        plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                                        p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                                        hold on;
                                        p.EdgeColor=Grey_lines;
                                       %plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                                            Site_no=(Extract{id,7});

                                         C=find(Cluster_members==Site_no);
                                         if isempty(C)
                                             TXT='';
                                         else
                                             TXT=Site_Channel_map{Site_no};
                                         end

                                       text((Extract{id,2}+1),-(Extract{id,3}+1),TXT,"FontSize",fontsize) 
                                                            
                                     end
                                  
                                    axis off
                                    axis equal
                                        
            %                         xl = xlim; % getting the X and Y limits to put in for later plots
            %                         yl= ylim;
                        
                                    hold on
            
            
            end
            
            
            
            function  Color_voronoi(example,color)
            global Extract
                for Site=1:1:size(example,1)
                    if example(Site) ==0
                        %do nothin
                    else
            
                        [coll locc]= find(cell2mat(Extract(:,7))==Site);
                            
                        if isempty(coll)
                                    for jj=1:1:size(Extract(:,9),1)
                                                Test_set=cell2mat(Extract(jj,9));
                                                if sum(Test_set==Site) > 0
                                                    Row=jj;
                                                end
                                              
                                    end
                        else
                            Row=coll;
                        end
            
                                               taco=Extract{Row, 6}(:,1);
                                               baco=Extract{Row, 6}(:,2);
                                
                                                p=fill(taco ,-baco,color,FaceAlpha=0.1);
                                                p.EdgeColor=color;
                    end
               end
            
            end
            
             function  Color_voronoi_Opac(example,color,Opac)
            global Extract
                for Site=1:1:size(example,1)
                    if example(Site) ==0
                        %do nothin
                    else
            
                        [coll locc]= find(cell2mat(Extract(:,7))==Site);
                            
                        if isempty(coll)
                                    for jj=1:1:size(Extract(:,9),1)
                                                Test_set=cell2mat(Extract(jj,9));
                                                if sum(Test_set==Site) > 0
                                                    Row=jj;
                                                end
                                              
                                    end
                        else
                            Row=coll;
                        end
            
                                               taco=Extract{Row, 6}(:,1);
                                               baco=Extract{Row, 6}(:,2);
                                
                                                p=fill(taco ,-baco,color,FaceAlpha=Opac);
                                                p.EdgeColor=color;
                    end
               end
            
            end
            
            
            
            function color = color_pick(number)
                % Define the colors in RGB format
            colors = [
                0.0, 0.0, 1.0;    % Blue
                1.0, 0.0, 0.0;    % Red
                0.0, 1.0, 0.0;    % Green
                1.0, 1.0, 0.0;    % Yellow
                0.0, 1.0, 1.0;    % Cyan
                1.0, 0.0, 1.0;    % Magenta
                0.5, 0.5, 0.5;    % Gray
                1.0, 0.5, 0.0;    % Orange
                0.5, 0.0, 0.5;    % Purple
                0.0, 0.5, 0.5;    % Teal
                0.5, 0.0, 0.0;    % Dark Red
                0.0, 0.5, 0.0;    % Dark Green
                0.0, 0.0, 0.5;    % Dark Blue
                0.75, 0.75, 0.0;  % Olive
                0.75, 0.0, 0.75;  % Violet
                0.0, 0.75, 0.75;  % Aqua
                0.5, 0.25, 0.0;   % Brown
                0.25, 0.5, 0.25;  % Light Green
                0.25, 0.25, 0.75; % Steel Blue
                0.75, 0.75, 0.75; % Light Gray
                0.0, 0.0, 0.0;    % Black
                1.0, 0.6, 0.0;    % Amber
                0.6, 1.0, 0.2;    % Lime Green
                0.6, 0.2, 1.0;    % Lavender
                1.0, 0.8, 0.6;    % Peach
                0.6, 0.8, 1.0;    % Sky Blue
                1.0, 0.2, 0.6;    % Rose
                0.4, 0.4, 0.8;    % Periwinkle
                0.8, 0.4, 0.4;    % Coral
                0.2, 0.6, 0.8;    % Turquoise
                0.8, 0.6, 0.2;    % Mustard
            ];

                % Ensure the number is within the valid range
                if number < 1 || number > size(colors, 1)
                    error('Input number must be between 1 and %d.', size(colors, 1));
                end
                
                % Output the corresponding color
                color = colors(number, :);
            end
    end %static
 end %Master function







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


