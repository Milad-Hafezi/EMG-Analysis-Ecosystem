
%% run the plotter
Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform

clearvars -except Meta_Data Unblank_raw mdl Extract
global Meta_Data Unblank_raw Extract mdl

close all
tic
mkdir MiladClustering\
mkdir Topoclustering_output\


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Plotters=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load the Meta Data file
if exist('Meta_Data')
else

[Partfile,Partpath] = uigetfile('*.mat','load the MetaData files');
                Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Channel_name);

end
%% Find the amplitude levels and sites
%since it should be similar for all channels, I only considered the first
%channel

       Test  = Meta_Data{1, 1}.Guide_Matrix;
             Unique_Sites= unique(Test(:,1));
             Unique_Amps=unique(Test(:,6));
             Channels=size(Meta_Data,2);
             Case=Meta_Data{1, 1}.Case;
             disp('Detected Amplitudes:  ')
             disp(Unique_Amps)
             Snips_fs=Meta_Data{1, 1}.Snips_fs; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for mama=1:1:16   % [3 4 2 8 9]
Chan=mama;

Amp=320;

Filename=sprintf('MiladClustering\\MiladCLusters(%s)_Ch%g.pdf',Case,Chan);
%% finding all mean for all the sites for each muscle
for j=1:1:size(unique(Meta_Data{1, 1}.Guide_Matrix(:,1)),1)
StructB=MiladFuncV1.Get_Train(Chan,j,Amp);
Mean_site_waveform(j,:)= trimmean(StructB.Trains,40,1); %mean(StructB.Trains); 
Mean_Info_site     (j,:)=StructB.Train_info(1,13);

end
Mean_struct.Trains=Mean_site_waveform;
Mean_struct.Trace='mean for all sites Ch8 256';
Mean_struct.Thresh=Mean_Info_site;
%% determining the testing set to perfom different DRMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing_Channel_means=Mean_struct.Trains;

if Unblank_raw ==1 
%         Testing_Channel_means=normalize(Mean_struct.Trains,'range',[0,1]);
%         Testing_Channel_means=Testing_Channel_means- mean(Testing_Channel_means,2); %offset

        Testing_Channel_means=Mean_struct.Trains-mean(Mean_struct.Trains,2);
else
        Testing_Channel_means=Mean_struct.Trains;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Blanking
if Unblank_raw ==1 
    Testing_Channel_means =MiladFuncV1.Blank_it_Train(Testing_Channel_means);
end

%% Filtration
%    Testing_Channel_means=MiladFuncV1.filtera(Testing_Channel_means,200);

%% Chop only post stim

%   Testing_Channel_means=Testing_Channel_means(:,1526:end);  % getting only post stim

%% Horizontal normalization
%      Testing_Channel_means=normalize(Testing_Channel_means,2,'range',[0,1]);

  %% normalization
%   if Unblank_raw ==1 
%       Testing_Channel_means=normalize(Testing_Channel_means,1,'range',[-0.5,0.5]);
%   else
%       Testing_Channel_means=normalize(Testing_Channel_means,1,'range',[0,1]);
%   end
% %% Take out the negative sites
%    Testing_Channel_means=Testing_Channel_means.*Mean_Info_site;

  
% Testing_set=Testing_set.*Mean_Info_site;

% Testing_set = filtera(Testing_set,30);
% [Testing_set,Cen,Scale] = normalize(abs(Testing_set), 2, "range",[0 1]);

% [Testing_set,Cen,Scale] = normalize(Testing_set, 2, "range",[-1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chopping to three windows
Testing_set=                Testing_Channel_means;
Testing_set_background=     Testing_set(:,1:1526);
Testing_set_stim=           Testing_set(:,1527:3053);
Testing_set_Post_stim=      Testing_set(:,3054:end);
Testing_set_POST_ONSET=[Testing_set_stim, Testing_set_Post_stim ];





%% Wavelet Signal Decomposition

% EMG_Signal=Testing_set_POST_ONSET;
% 
% clear Decomposed_Signal_High Decomposed_Signal_Mid Decomposed_Signal_Low
% for i=1:1:size(EMG_Signal,1)
% 
%     Signal=EMG_Signal(i,:);
%     levelForReconstruction_High= [false, false,  true, false, false, false, false];
%     levelForReconstruction_Mid=  [false, false, false, false, false, true,  false];
%     levelForReconstruction_Low=  [false, false, false, false, false, false,  true];
%     % Perform the decomposition using modwt
%     wt = modwt(Signal, 'sym2', 6);
%     % Construct MRA matrix using modwtmra
%     mra = modwtmra(wt, 'sym2');
%     % Sum along selected multiresolution signals
%     Decomposed_Signal_High(i,:) = sum(mra(levelForReconstruction_High,:),1);
%     Decomposed_Signal_Mid (i,:) = sum(mra(levelForReconstruction_Mid,:) ,1);
%     Decomposed_Signal_Low (i,:) = sum(mra(levelForReconstruction_Low,:) ,1);
% end


% Reassignement
% Testing_set_POST_ONSET=Decomposed_Signal_Low;

%% Dynamic Time Warping Distance

% 
% maxWarpingWindow = 150;

% Compute DTW distance matrix
% DWT_Original =      MiladFuncV1.computeDTWMatrix(Testing_set, maxWarpingWindow);
% DWT_Post_Onset=     MiladFuncV1.computeDTWMatrix(Testing_set_POST_ONSET, maxWarpingWindow);
% 
% DWT_High =          MiladFuncV1.computeDTWMatrix(Decomposed_Signal_High, maxWarpingWindow);
% DWT_Mid =           MiladFuncV1.computeDTWMatrix(Decomposed_Signal_Mid , maxWarpingWindow);
% DWT_Low =           MiladFuncV1.computeDTWMatrix(Decomposed_Signal_Low , maxWarpingWindow);
% 

% beep
% 
% 
% 
% %% Plotting & deduced measure for original
% 
% Deduced_measure= sum(DWT_Post_Onset,1);

%%
Logic_verdict=(Mean_Info_site>0);
Target=Testing_set_POST_ONSET;

% Testing_set_POST_ONSET=normalize(Testing_set_POST_ONSET,2,"range",[0 1]);


% Target=MiladFuncV1.filtera(Testing_set_POST_ONSET,45);

Target=Target.*Logic_verdict;

[activeIndices, activeSignals] = MiladFuncV1.separateActiveSignals(Target);

%% Miladian clustering

% Inputs:
% Testing_set_POST_ONSET: NxM matrix of EMG waveforms (N waveforms, M samples each)
% maxDepth: maximum depth for recursion
% minClusterSize: minimum cluster size to allow further splitting
% maxWarpingWindow: parameter for your DTW calculation

maxDepth = 5;  %8;        % example
minClusterSize = 3; % example
maxWarpingWindow = 150; % example


%  activeSignals=normalize(activeSignals,2,"range",[0 1]);
% Example usage after hierarchical clustering:
%  finalClusters = MiladFuncV1.hierarchical_linkage_with_recompute(Target, maxDepth, minClusterSize, maxWarpingWindow);

 [finalClusters, progressReport] = MiladFuncV1.hierarchical_linkage_with_recompute(activeSignals, maxDepth, minClusterSize, maxWarpingWindow);

 % disp(progressReport)

 Revised_report=MiladFuncV1.plotClusteringTree(progressReport,activeIndices);
 if Plotters==1
   exportgraphics(gca,Filename,'Append',true)
     pause(1);
 end
 
  

PreClusters=ones(size(Target,1),1);
Uped=finalClusters+1;
PreClusters(activeIndices)=Uped;


clear MEMBERCELL_UnCleaned
for i = unique(PreClusters)'
    clusterID = i;
    clusterMembers = find(PreClusters == clusterID); % Get EMG indices in this cluster
    MEMBERCELL_UnCleaned{i,1}=Testing_Channel_means(clusterMembers,:); 
end

%  Clusters= finalClusters;

Master_Cluster{Amp,1,Chan}=PreClusters; % level one clusters

%% Now reassign small clusters:

MemberLimit=2; % goes and grabs any cluster with less than "MemberLimit" number of elements
maxWindow=150;
% max_Delay=50;
Reduced_cluster = MiladFuncV1.reassignSmallClusters(finalClusters, activeSignals, MemberLimit,maxWindow);

Sweep_Clusters=ones(size(Target,1),1);
Uped=Reduced_cluster+1;
Sweep_Clusters(activeIndices)=Uped;

Master_Cluster{Amp,2,Chan}=Sweep_Clusters; % level two clusters

%% finding latency

 baselineCols = 1:1526;

 
    % 2) Compute baseline stats
    [meanVals, stdVals] = MiladFuncV1.getBaselineStats(Testing_set, baselineCols);
    disp('Baseline means:');
    disp(meanVals);
    disp('Baseline std devs:');
    disp(stdVals);

    % 3) Double-threshold method
    lowFactor  = 2;
    highFactor = 5;
    minPersistSamples = 10;
    doubleThreshOnsets = MiladFuncV1.doubleThresholdOnsetAll(...
        Testing_set, lowFactor, highFactor, baselineCols, minPersistSamples);
    
%     nnz(Logic_verdict) - nnz(Logic_verdict'.*(~isnan(doubleThreshOnsets)))
%     MiladFuncV1.plot_latency(Testing_set,doubleThreshOnsets)


 

    % 6) Adaptive threshold (All correct)
    factor     = 2;
    minDuration = 15;
    adapOnsets = MiladFuncV1.adaptiveThresholdAll(Testing_set, ...
        baselineCols, factor, minDuration);
   
%     nnz(Logic_verdict) - nnz(Logic_verdict'.*(~isnan(adapOnsets)))
%     MiladFuncV1.plot_latency(Testing_set,adapOnsets)



    % 8) Wavelet-based onset (doing really good but missing many signals (Nan))

    waveletName = 'db4';
    level       = 4;
    waveFactor  = 2;
    minDur      = 10;
    waveOnsets  = MiladFuncV1.waveletOnsetAll(Testing_set, baselineCols, ...
        waveletName, level, waveFactor, minDur);

%    nnz(Logic_verdict) - nnz(Logic_verdict'.*(~isnan(waveOnsets)))

%      MiladFuncV1.plot_latency(Testing_set,waveOnsets)


%% Now find the clusters that are large latency
 % Define a factor for deciding “significantly different”
    % e.g., if average lag > factor * (cluster avg-lag) => outlier

    

% maxLag=500;
% updatedCluster = MiladFuncV1.reassignByAvgLag(activeSignals, Reduced_cluster, maxLag);
% 
% aa= Reduced_cluster-updatedCluster
% 
% Clusters=ones(size(Target,1),1);
% Uped=updatedCluster+1;
% Clusters(activeIndices)=Uped;

%% In case there is a tiny signal, I do this using the normalization


maxLag_Norm=500;
Latency_vector_active=doubleThreshOnsets (activeIndices); % adapOnsets   % waveOnsets %doubleThreshOnsets

normed_Signal=normalize(activeSignals,2,'range',[0 1]);
Normed_updatedCluster_initial = MiladFuncV1.reassignByAvgLag(normed_Signal, Reduced_cluster, Latency_vector_active, maxLag_Norm);

        % redo the small clusters
        Reduced_cluster2 = MiladFuncV1.reassignSmallClusters(Normed_updatedCluster_initial, activeSignals, MemberLimit,maxWindow);
        Normed_updatedCluster = MiladFuncV1.reassignByAvgLag(normed_Signal, Reduced_cluster2, Latency_vector_active, maxLag_Norm);





% a= Reduced_cluster-Normed_updatedCluster

Clusters=ones(size(Target,1),1);
Uped=Normed_updatedCluster+1;
Clusters(activeIndices)=Uped;








%% mapping the clusters
% Get the "stable" unique values and an index vector:
[uniqueVals, ~, idxMap] = unique(Clusters, 'stable');

% idxMap now holds the integer labels for each element of v.
% By default, uniqueVals = [1   2   5   11],
% and idxMap =           [1   1   2   2   3    4   4].

% If you just want consecutive integers:
Clusters = idxMap;
%%



Master_Cluster{Amp,3,Chan}=Clusters; % Last level clusters
Master_Latency{Amp,Chan}=doubleThreshOnsets; % level two clusters
%% ploting the cleaned clusters
if Plotters==1
        clear MEMBERCELL_Cleaned
        for i = unique(Clusters)'
            clusterID = i;
            clusterMembers = find(Clusters == clusterID); % Get EMG indices in this cluster
            MEMBERCELL_Cleaned{i,1}=Testing_Channel_means(clusterMembers,:);
        
        %     subplot(length(uniqueClusters), 1, i); % Create a subplot for this cluster
            F=figure;
            for j = 1:length(clusterMembers)
                subplot(3,1,1)
                    if 1  % Mean_Info_site(clusterMembers(j))>0
                        plot(Testing_Channel_means(clusterMembers(j), :),"DisplayName",string(clusterMembers(j))); % Plot each EMG in this cluster
                        title(sprintf('Original Cluster %s' , num2str(clusterID)))
                         hold on;
                         legend 
                    end
                subplot(3,1,2)
                        plot(Target(clusterMembers(j), :),"DisplayName",string(clusterMembers(j))); % Plot each EMG in this cluster
                        title(sprintf('Target Cluster %s' , num2str(clusterID)))
                         hold on;
                         legend 
                subplot(3,1,3)  
                        plot(normalize(Target(clusterMembers(j), :),2,'range',[0 1]),"DisplayName",string(clusterMembers(j))); % Plot each EMG in this cluster
                        title(sprintf('normalized-filter 45 Cluster %s' , num2str(clusterID)))
                         hold on;
                         legend 
        
           
                         
            end
            hold off;
            
            sgtitle(sprintf('(%s) Win1:%g Win2:%g #Clust:%g CH:%g AMP:%g',Case,maxWarpingWindow,maxWindow,size(unique(Clusters), 1),Chan,Amp));
            xlabel('Time (Samples)');
            ylabel('EMG Signal');
            grid on;
             exportgraphics(F,Filename,'Append',true)
             pause(1);
        
        
        end



        %% ploting clean trends
        
        FF=figure;
        FF.WindowState = 'maximized'; % Set the figure to full screen
        
        for J=1:1:size(MEMBERCELL_Cleaned,1)
                 Trend=median(MEMBERCELL_Cleaned{J,1},1);
                 Trend=normalize(Trend, 2, 'range',[0 1]);
        
        %          subplot(size(MEMBERCELL_Cleaned,1),1,J);
                 plot(Trend+J,"DisplayName",num2str(J));
                 hold on
        %          title(sprintf('Cluster# %g Trend cleaned',J))
        end
        sgtitle(sprintf('Cleaned_and_Latencied (%s) Win1:%g Win2:%g #Clust:%g CH:%g AMP:%g',Case,maxWarpingWindow,maxWindow,size(unique(Clusters), 1),Chan,Amp));
        legend
        exportgraphics(FF,Filename,'Append',true)
        pause(1);
        
        
        %% ploting UNclean trends
        FFF=figure;
        FFF.WindowState = 'maximized'; % Set the figure to full screen
        
        
        for J=1:1:size(MEMBERCELL_UnCleaned,1)
                 Trend=median(MEMBERCELL_UnCleaned{J,1},1);
                 Trend=normalize(Trend, 2, 'range',[0 1]);
        
        %          subplot(size(MEMBERCELL_Cleaned,1),1,J);
                 plot(Trend+J,"DisplayName",num2str(J));
                 hold on
        %          title(sprintf('Cluster# %g Trend cleaned',J))
        end
        sgtitle(sprintf('UNCLEANED (%s) Win1:%g Win2:%g #Clust:%g CH:%g AMP:%g',Case,maxWarpingWindow,maxWindow,size(unique(Clusters), 1),Chan,Amp));
        legend
        exportgraphics(FFF,Filename,'Append',true)
        pause(1);



end



%% plotting clusters
% sgtitle('EMG Signals by Cluster');

for Tac=1:1:2
if Tac==1
    Clusters_of_interest=PreClusters;
else
    Clusters_of_interest=Clusters;
end


numClusters=length(unique(Clusters_of_interest));



cluster_groups = cell(numClusters, 1);
for idx = 1:length(Clusters_of_interest)
    cluster_groups{Clusters_of_interest(idx)} = [cluster_groups{Clusters_of_interest(idx)}, idx];
end




% then I need to have a code or function to plot and color the clusters

Verd=(Mean_Info_site>0);

Verdictor(:,Chan)=(Mean_Info_site>0 & Mean_Info_site<=Amp);

if Plotters==1
                h = zeros(1, size(cluster_groups, 1));
                MiladFuncV1.plot_voronoi() % plotting the voronois
                
                hold on
                
                for clst=1:1:size(cluster_groups,1)
                    color = MiladFuncV1.color_pick(clst);
                
                    Preload=zeros(size(Testing_Channel_means,1),1);
                    Preload(cluster_groups{clst})=1;
                    Final=Preload.*Verd;
                
                    MiladFuncV1.Color_voronoi(Final,color) % be careful to only feed in the column for each cluster
                    color_used(clst,:)=color;
                    h(clst) = plot(0, 0, 's', 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k'); % Invisible marker to create legend entry
                end
                
                % Generate the legend dynamically with the correct colors and labels
                legend(h, arrayfun(@(x) sprintf('Cluster %d', x), 1:size(cluster_groups, 1), 'UniformOutput', false), 'Location', 'Best');
                
                
                note='';
                
                if Tac==1
                    
                    title(sprintf('UNCLEANED PreClusters \n Win:%g #Clust:%g CH:%g AMp:%g \n  %s',maxWarpingWindow,size(cluster_groups, 1),Chan,Amp,note))
                else
                   
                    title(sprintf('Cleaned_and_Latencied \n Win:%g #Clust:%g CH:%g AMp:%g \n  %s',maxWarpingWindow,size(cluster_groups, 1),Chan,Amp,note))
                end
                
                
                exportgraphics(gca,Filename,'Append',true)
                pause(1);
                % 
                beep
                
                
                
                close all
end


        % Building the muscle signituare waveforms:


            if Tac==1
                     for clst=1:1:size(cluster_groups,1)
                        
                        Cluster_memebers_Pre{clst,Chan}=Testing_Channel_means(cluster_groups{clst},:);
                        Clusters_Pre{clst,Chan}=cluster_groups{clst};
                     end
                            
                        
                 
  
           
            elseif Tac==2
                 


                       for clst=1:1:size(cluster_groups,1)
                        
                        Cluster_memebers_Refined{clst,Chan}=Testing_Channel_means(cluster_groups{clst},:);
                        Clusters_Refined{clst,Chan}=cluster_groups{clst};
            
                       end




            end
            

end




end % for mama
  

[Signature_Pre,Tracker_Pre,Combined_Cluster_ident_Pre,Combined_Verdict_Pre]=MiladFuncV1.properCluster(Cluster_memebers_Pre,Clusters_Pre,Verdictor);

[Signature_Refined,Tracker_Refined,Combined_Cluster_ident_Refined,Combined_Verdict_Refined]=MiladFuncV1.properCluster(Cluster_memebers_Refined,Clusters_Refined,Verdictor);



% save(sprintf('Topoclustering_output/Topo_clusters(%s)_%s_Raw%g.mat',Case,datestr(now, 'dd-mmm-yyyy'),Unblank_raw),'Master_Cluster','Master_Latency','Cluster_memebers_Pre','Cluster_memebers_Refined','Clusters_Refined','Clusters_Pre','Tracker_Pre', 'Combined_Cluster_ident_Pre','Tracker_Refined','Combined_Cluster_ident_Refined','Unblank_raw','Signature_Pre','Signature_Refined','Combined_Verdict_Pre','Combined_Verdict_Refined','Case','Snips_fs','-v7.3')

% First, create the struct that will hold everything
Topo_Muscle_Cluster.Master_Cluster              = Master_Cluster;
Topo_Muscle_Cluster.Master_Latency              = Master_Latency;
Topo_Muscle_Cluster.Cluster_memebers_Pre        = Cluster_memebers_Pre;
Topo_Muscle_Cluster.Cluster_memebers_Refined    = Cluster_memebers_Refined;
Topo_Muscle_Cluster.Clusters_Refined            = Clusters_Refined;
Topo_Muscle_Cluster.Clusters_Pre                = Clusters_Pre;
Topo_Muscle_Cluster.Tracker_Pre                 = Tracker_Pre;
Topo_Muscle_Cluster.Combined_Cluster_ident_Pre  = Combined_Cluster_ident_Pre;
Topo_Muscle_Cluster.Tracker_Refined             = Tracker_Refined;
Topo_Muscle_Cluster.Combined_Cluster_ident_Refined = Combined_Cluster_ident_Refined;
Topo_Muscle_Cluster.Unblank_raw                 = Unblank_raw;
Topo_Muscle_Cluster.Signature_Pre               = Signature_Pre;
Topo_Muscle_Cluster.Signature_Refined           = Signature_Refined;
Topo_Muscle_Cluster.Combined_Verdict_Pre        = Combined_Verdict_Pre;
Topo_Muscle_Cluster.Combined_Verdict_Refined    = Combined_Verdict_Refined;
Topo_Muscle_Cluster.Case                        = Case;
Topo_Muscle_Cluster.Snips_fs                    = Snips_fs;
Topo_Muscle_Cluster.Size                        = size(Testing_Channel_means);


% Now save this single struct with -v7.3
save(sprintf('Topoclustering_output/Topo_clusters(%s)_%s_Raw%g.mat', ...
    Case, datestr(now, 'dd-mmm-yyyy'), Unblank_raw), ...
    'Topo_Muscle_Cluster', '-v7.3');
