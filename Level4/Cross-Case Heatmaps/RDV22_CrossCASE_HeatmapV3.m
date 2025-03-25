

% clearvars -except Meta_Data Unblank_raw mdl Topo_Muscle_Cluster Extract mdl
clear all

global Meta_Data Unblank_raw  Extract mdl

close all
tic
mkdir CrossCASE\


% mkdir CrossMuscles_Topoclustering_output\




%% Batch upload


Case_set=["22-03", "22-04", "22-05", "22-06", "22-07"];

Super_Signatures=[];
Cases=[];
All_Number_of_Sites=[];
Super_Case_Tracker=[];
Super_Tracker=[];
Super_Logic_verdict=[];
Current_index=0;


for J=1:1:length(Case_set)

                     [Partfile,Partpath] = uigetfile('*.mat',sprintf('load Topo files for %s',Case_set(J)));
                    Load_Channel_name=sprintf('%s\\%s',Partpath,Partfile);
                    load(Load_Channel_name);

          
          fieldName = ['Case' num2str(J)];


         Superfile_complete.(fieldName)=Topo_Muscle_Cluster;  
        
         Superfile.(fieldName).Signatures=                 Topo_Muscle_Cluster.Signature_Refined;
         Superfile.(fieldName).Combined_Cluster_ident=     Topo_Muscle_Cluster.Combined_Cluster_ident_Refined;
         Superfile.(fieldName).Tracker =                   Topo_Muscle_Cluster.Tracker_Refined;
         Superfile.(fieldName).Logic_verdict=              (Topo_Muscle_Cluster.Combined_Verdict_Refined > 0);
         Superfile.(fieldName).Case =                      Topo_Muscle_Cluster.Case;
         Superfile.(fieldName).Snips_fs =                  Topo_Muscle_Cluster.Snips_fs;
         Superfile.(fieldName).Number_of_Sites =           Topo_Muscle_Cluster.Size (1);    
          Superfile.(fieldName).Unblank_raw =             Topo_Muscle_Cluster.Unblank_raw;    
         

        %building the mass matrix
        Counter=0;
        for i=1:1:size(Superfile.(fieldName).Signatures,2)
                Super_Signatures=[Super_Signatures;  Superfile.(fieldName).Signatures{i}];
                Counter=Counter+ size(Superfile.(fieldName).Signatures{i},1);
        end
        

        Cases{J,1}=Superfile.(fieldName).Case;
        

        for k=1:1:Counter
                
                 Super_Combined_Cluster_ident{Current_index+k,1} = Superfile.(fieldName).Combined_Cluster_ident{k};
        end
                Current_index=Current_index+Counter;

         
        
        Super_Case_Tracker=[Super_Case_Tracker ;(ones(Counter,1) * J)];
        Super_Tracker= [Super_Tracker; Superfile.(fieldName).Tracker];

        Super_Logic_verdict=[Super_Logic_verdict; Superfile.(fieldName).Logic_verdict];


        All_Number_of_Sites(J)=Superfile.(fieldName).Number_of_Sites;
        Raw_count(J)=Superfile.(fieldName).Unblank_raw;
        
end
%%  Defininf inputs
if sum(Raw_count) > 0
    Unblank_raw=1;
else
    Unblank_raw=0;
end


Filename=sprintf('CrossCASE\\CrossCASE_(%s)_Raw%g.pdf',datestr(now, 'dd-mmm-yyyy'),Unblank_raw);

%% Outputs so far

Cases;
All_Number_of_Sites;



Super_Signatures;
Super_Logic_verdict;
Super_Tracker;
Super_Combined_Cluster_ident; 



Super_Case_Tracker;





%% Reassigning the new variables to match old code
Logic_verdict=Super_Logic_verdict;
Tracker=Super_Tracker;
Combined_Cluster_ident=Super_Combined_Cluster_ident;

%% determining the testing set to perfom different DRMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing_Channel_means=Mean_struct.Trains;

if Unblank_raw ==1 
%         Testing_Channel_means=normalize(Mean_struct.Trains,'range',[0,1]);
%         Testing_Channel_means=Testing_Channel_means- mean(Testing_Channel_means,2); %offset

        Testing_Channel_means=Super_Signatures-mean(Super_Signatures,2);
else
        Testing_Channel_means=Super_Signatures;
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


%% Define Target
Target= Testing_set;%Testing_set;  %Decomposed_Signal_Low;

Target=normalize(Target,2,"range",[0 1]);

% Target=MiladFuncV1.filtera(Testing_set_POST_ONSET,45);

Target=Target.*Logic_verdict;
[activeIndices, activeSignals] = MiladFuncV1.separateActiveSignals(Target);


Reduced_Tracker=Tracker(activeIndices,:);
Reduced_Combined_Cluster_ident=Combined_Cluster_ident(activeIndices,:);
Reduced_Super_Case_Tracker=Super_Case_Tracker(activeIndices,:);


%% PCs disection


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

[Low_order_approximation,SuperLeft,SuperRight]=MiladFuncV1.Milad_NMF(NMF_Target,Factors,algorithm);






%% Score ratio for each factor based on Case
clear Factor_Score Case_Ratio
for i=1:1:size(SuperRight,1)
    Factor_Score(i,1)=sum(SuperLeft(:,i));

    for J=unique(Reduced_Super_Case_Tracker)'
          Case_Ratio(i,J)=   sum(SuperLeft(Reduced_Super_Case_Tracker==J,i))/sum(SuperLeft(:,i)) ;
    end


end



%% Strategy for vornoi plots

% one figure with the trend

% A speperate figure for all the heatmaps.
            

%% Breaking the data for each case
for Cas=unique(Reduced_Super_Case_Tracker)'
        clear Locator
        Locator= Reduced_Super_Case_Tracker ==Cas;

            Left=SuperLeft(Locator>0,:);
            Right=SuperRight ;
            Number_of_Sites=All_Number_of_Sites(Cas);
            
            Sub_Reduced_Tracker=Reduced_Tracker(Locator,:);

            g=find(Locator>0);
            clear Sub_Reduced_Combined_Cluster_ident
            for S=1:1:length(g)
               Sub_Reduced_Combined_Cluster_ident{S,1}=Reduced_Combined_Cluster_ident{g(S)};
            end

            
            
            
            
            %% Channel_waveform score (Per Case)
            
            Norm_Left=normalize(Left,1,'range',[0 1]);
            % Norm_Left=Left;
            clear Channel_Score
            for H=1:1:size(Norm_Left,2)
                   Scores =Norm_Left(:,H);
                    for i=1:1:16
                        Channel_Score(i,H)=sum(Scores(Sub_Reduced_Tracker(:,1)==i));
                    end
            end
            
            
            %% Site-Waveform Score  (Per Case)
            
            Norm_Left=normalize(Left,1,'range',[0 1]);
            
            clear Site_Score
            for H=1:1:size(Norm_Left,2)
                   Scores =Norm_Left(:,H);
                    for i=1:1:Number_of_Sites
                            clear mapper
                            for j=1:1:length(Sub_Reduced_Combined_Cluster_ident)
                                 mapper(j)=sum(Sub_Reduced_Combined_Cluster_ident{j}==i);
                            end
            
                            Site_Score(i,H)=sum(Scores(mapper>0));
                    end
            end
            
            
            
           %% warning
           DIA=msgbox(sprintf('Make sure to select the voronoi for case %s',Cases{Cas}),'WARNING!');
           uiwait(DIA)
                       
            %% plotting clusters in voronois
            % sgtitle('EMG Signals by Cluster');
            
            %
            Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform_CrossCase
            for clst=1:1:size(Site_Score, 2)
                Figgs=figure;
                Figgs.WindowState = 'maximized'; % Set the figure to full screen
            
                
            %     Base_color = MiladFuncV1.color_pick(clst);
                fontsize=2;
                          
            
                MiladFuncV1.Subplot_voronoi([],[],fontsize); % plotting the voronois
                hold on
            
            
                Norm_scoerd=Case_Ratio(clst,Cas)*normalize(Site_Score(:,clst),1,'range',[0 1]);
                
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
                     
                
                 
                
                  sgtitle(sprintf('%s Hotsposts for PC %g',Cases{Cas},clst))
                  exportgraphics(Figgs,Filename,'Append',true)
            end
            
            
            close all



end


%%






% 
% 
% 
% save(sprintf('Topoclustering_output/Topo_clusters(%s)_%s_Raw%g.mat',Case,datestr(now, 'dd-mmm-yyyy'),Unblank_raw),'Master_Cluster_Cross_muscles','Master_Latency','Cluster_memebers_Pre','Cluster_memebers_Refined','Clusters_Refined','Clusters_Pre','Unblank_raw','-v7.3')
% 
% 

