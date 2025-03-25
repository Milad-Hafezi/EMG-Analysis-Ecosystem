

Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform
clearvars -except Meta_Data Unblank_raw mdl Topo_Muscle_Cluster Extract mdl
global Meta_Data Unblank_raw  Extract mdl

close all
tic
mkdir CrossMuscles\
% mkdir CrossMuscles_Topoclustering_output\
%% load the Meta Data file
if exist('Meta_Data')
else
    [Partfile,Partpath] = uigetfile('*.mat','load the MetaData files');
                    Load_Channel_name=sprintf('%s\\%s',Partpath,Partfile);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('Topo_Muscle_Cluster')
else

    [Partfile,Partpath] = uigetfile('*.mat','load the TopoCluster files');
                    Load_Channel_name=sprintf('%s\\%s',Partpath,Partfile);
                    load(Load_Channel_name);
end

%% Extracting the signautures for each muscle (this could be mean, median or anything else)

Signatures=                 Topo_Muscle_Cluster.Signature_Refined;
Combined_Cluster_ident=     Topo_Muscle_Cluster.Combined_Cluster_ident_Refined;
Tracker =                   Topo_Muscle_Cluster.Tracker_Refined;
Logic_verdict=              (Topo_Muscle_Cluster.Combined_Verdict_Refined > 0);


Number_of_Sites=Topo_Muscle_Cluster.Size (1);


%% building the combined Matrix

% Initialize combined matrix and a tracker matrix
Combined_Muscles_EMGs = [];
for i = 1:length(Signatures)
    currentBlock = Signatures{i}; 
    Combined_Muscles_EMGs = [Combined_Muscles_EMGs; currentBlock];
    
end


%% outputs so far


Combined_Muscles_EMGs;
Tracker;
Combined_Cluster_ident; 
Logic_verdict; %Combined_Verdict;



 Filename=sprintf('CrossMuscles\\CrossMuscle_(%s).pdf',Case);



%% determining the testing set to perfom different DRMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing_Channel_means=Mean_struct.Trains;

if Unblank_raw ==1 
%         Testing_Channel_means=normalize(Mean_struct.Trains,'range',[0,1]);
%         Testing_Channel_means=Testing_Channel_means- mean(Testing_Channel_means,2); %offset

        Testing_Channel_means=Combined_Muscles_EMGs-mean(Combined_Muscles_EMGs,2);
else
        Testing_Channel_means=Combined_Muscles_EMGs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Blanking
if Unblank_raw ==1 
    Testing_Channel_means =MiladFuncV1.Blank_it_Train(Testing_Channel_means);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Chopping to three windows
Testing_set=                Testing_Channel_means;
Testing_set_background=     Testing_set(:,1:1526);
Testing_set_stim=           Testing_set(:,1527:3053);
Testing_set_Post_stim=      Testing_set(:,3054:end);
Testing_set_POST_ONSET=[Testing_set_stim, Testing_set_Post_stim ];

%% Wavelet Signal Decomposition

% EMG_Signal=Testing_set;
% 
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
% 
%% filtration
% Filtered=MiladFuncV1.filtera(Testing_set,45);
% 


%% Define Target
Target= Testing_set; %Testing_set;  %Decomposed_Signal_Low;

Target=normalize(Target,2,"range",[0 1]);

% Target=MiladFuncV1.filtera(Testing_set_POST_ONSET,45);

Target=Target.*Logic_verdict;
[activeIndices, activeSignals] = MiladFuncV1.separateActiveSignals(Target);


Reduced_Tracker=Tracker(activeIndices,:);
Reduced_Combined_Cluster_ident=Combined_Cluster_ident(activeIndices,:);




%% PCs disection


% NMF_Target=Testing_set;  % activeSignals  %Target
% Reduced_Tracker=Tracker;
% Reduced_Combined_Cluster_ident=Combined_Cluster_ident;



NMF_Target=activeSignals;  % activeSignals  %Target

% NMF_Target=NMF_Target(:,1527:3053);

NMF_Target=NMF_Target(:,1527:end);

maxFactors=20;
optimalFactors=MiladFuncV1.nnmf_Nfactor_error(NMF_Target,maxFactors)

 Factors=optimalFactors;

%  if Factors>10
%  Factors=10;
%  end

algorithm="mult";   % 'mult'; % "als";

[Low_order_approximation,Left,Right]=MiladFuncV1.Milad_NMF(NMF_Target,Factors,algorithm);


%% Channel_waveform score

Norm_Left=normalize(Left,1,'range',[0 1]);
% Norm_Left=Left;
clear Channel_Score
for H=1:1:size(Norm_Left,2)
       Scores =Norm_Left(:,H);
        for i=1:1:16
            Channel_Score(i,H)=sum(Scores(Reduced_Tracker(:,1)==i));
        end
end


%% Site-Waveform Score




Norm_Left=normalize(Left,1,'range',[0 1]);

clear Site_Score
for H=1:1:size(Norm_Left,2)
       Scores =Norm_Left(:,H);
        for i=1:1:Number_of_Sites
                clear mapper
                for j=1:1:length(Reduced_Combined_Cluster_ident)
                     mapper(j)=sum(Reduced_Combined_Cluster_ident{j}==i);
                end

                Site_Score(i,H)=sum(Scores(mapper>0));
        end
end







%% plotting clusters in voronois
% sgtitle('EMG Signals by Cluster');




% then I need to have a code or function to plot and color the clusters

% Verd=(Mean_Info_site>0);
% h = zeros(1, size(Site_Score, 2));
% MiladFuncV1.plot_voronoi() % plotting the voronois

% hold on
%
for clst=1:1:size(Site_Score, 2)
    Figgs=figure;
    Figgs.WindowState = 'maximized'; % Set the figure to full screen

    
%     Base_color = MiladFuncV1.color_pick(clst);
    fontsize=2;
    subplot(2,1,1)
    MiladFuncV1.Subplot_voronoi([],[],fontsize); % plotting the voronois
    hold on


    Norm_scoerd=normalize(Site_Score(:,clst),1,'range',[0 1]);

    for J=1:1:Number_of_Sites
           Preload=zeros(Number_of_Sites,1);
           Preload(J,:)=1;
           Final=Preload; %.*Verd;
           color= [0,1-Norm_scoerd(J),0]; % Green
%            color= [Norm_scoerd(J),0,0]; % Red

           FaceAlpha=1;   % 0.1;

           MiladFuncV1.Color_voronoi_Opac(Final,color,FaceAlpha) % be careful to only feed in the column for each cluster

    end


  subplot(2,1,2)

     plot(Right(clst,:))%,"DisplayName",string(cluster_groups{clst}))
         
    
     
    
      sgtitle(sprintf('%s Hotsposts for PC %g',Case,clst))
      exportgraphics(Figgs,Filename,'Append',true)
end





close all
%%






% 
% 
% 
% save(sprintf('Topoclustering_output/Topo_clusters(%s)_%s_Raw%g.mat',Case,datestr(now, 'dd-mmm-yyyy'),Unblank_raw),'Master_Cluster_Cross_muscles','Master_Latency','Cluster_memebers_Pre','Cluster_memebers_Refined','Clusters_Refined','Clusters_Pre','Unblank_raw','-v7.3')
% 
% 

