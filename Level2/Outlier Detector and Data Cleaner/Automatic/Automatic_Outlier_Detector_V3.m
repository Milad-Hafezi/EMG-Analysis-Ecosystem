
% 

tic

clear all
close all
import mlreportgen.report.*
import mlreportgen.dom.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PDF_logic=0; % Tunrs the PDF report generation On (1) or off (anything but 1)
% More info on PDF Reports
% https://www.mathworks.com/help/rptgen/ug/center_figure_snapshot_on_a_page.html
% https://www.mathworks.com/help/rptgen/ug/mlreportgen.dom.pdfpagelayout-class.html


Descriptor=''


% Blanking


%% load the Meta Data file
[Partfile,Partpath] = uigetfile('*.mat','load the MetaData files');
                Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Channel_name);

Version_file=Partfile;

%% Find the amplitude levels and sites
%since it should be similar for all channels, I only considered the first
%channel

       test  = Meta_Data{1, 1}.Guide_Matrix;
             Unique_Sites= unique(test(:,1));
             Unique_Amps=unique(test(:,6));
             Channels=size(Meta_Data,2);
             Number_of_Tirals=size(test,1);
             disp('Detected Amplitudes:  ')
             disp(Unique_Amps)

%% Making directories
mkdir AutoCleanedData

Selected_Train_Matrix=nan(Number_of_Tirals,Channels);

%% Loop in channels
for Selected_Channel=1:1:Channels
          % Assigne the proper Meta Data for the channel
          clear Guide_trains Guide_Matrix Muscle_name Case Channel_Number Snips_fs
          Guide_trains  = Meta_Data{1, Selected_Channel}.Guide_trains ; 
          Guide_Matrix  = Meta_Data{1, Selected_Channel}.Guide_Matrix;
          Muscle_name   = Meta_Data{1, Selected_Channel}.Muscle_name;
                  Case  = Meta_Data{1, Selected_Channel}.Case;
        Channel_Number  = Meta_Data{1, Selected_Channel}.Channel_Number;
              Snips_fs  = Meta_Data{1, Selected_Channel}.Snips_fs;
             
    
          % PDF Report setup
                if PDF_logic==1
                Rep = Report(sprintf('%s_Ch%g_%s',Case,Channel_Number,Muscle_name),'pdf');
                open(Rep);
                % setting up the report layout
                
                if strcmpi(Rep.Type,"pdf")
                    pageLayoutObj = PDFPageLayout;
                else
                    pageLayoutObj = DOCXPageLayout;
                end
                %Specify the page orientation, height, and width.
                pageLayoutObj.PageSize.Orientation = "portrait"; %landscape   %portrait
                pageLayoutObj.PageSize.Height = "11 in";
                pageLayoutObj.PageSize.Width = "8.5 in";
                
                %Specify the page margins.
                pageLayoutObj.PageMargins.Top = "0 in";
                pageLayoutObj.PageMargins.Bottom = "0 in";
                pageLayoutObj.PageMargins.Left = "0 in";
                pageLayoutObj.PageMargins.Right = "0 in";
                pageLayoutObj.PageMargins.Header = "0 in";
                pageLayoutObj.PageMargins.Footer = "0 in";
                
                %Add the page layout object to the report.
                add(Rep,pageLayoutObj);
                end

        %% loop in sites
        for Amp_interest=flip(sort(Unique_Amps))' % remember that you need to have a horizontal vecotr to pull this off
                        
                     [row,col, val] =find (Guide_Matrix(:,6)==Amp_interest);
                                      
                     clear Sub_matrix_Trains Sub_matrix_info
                     Sub_matrix_info= Guide_Matrix(row,:);
                     Sub_matrix_Trains=Guide_trains(row,:);

                            %% Amplitude loop
                            for     god=Unique_Sites'     % remember that you need to have a horizontal vecotr to pull this off
                                        clear Trains   Train_info
                                            % finding the correct sites
                                               [ro,co,v]=find (Sub_matrix_info(:,1)==god);
                                                    clear Bypass
                                                    if isempty(ro)
                                                        Bypass=1;

                                                    else
                                                        Train_info=Sub_matrix_info(ro,:);
                                                        Trains=Sub_matrix_Trains(ro,:);
                                                    end
%% Here is the core code
if exist('Bypass')
    disp('Amp not found for this site')
else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sectioning 
Precapture_index=Meta_Data{1, 1}.Guide_Matrix(1,9);

    


%% Blanking
Blanking_control=0;
Pulse_Period=4.99712;   %ms     0.4897  0.00499712
Blanking_window=3.8; % ms
Blanking_Shift=0; % a number from 0 to 100 (percentage) this can shift all the blanking windows by a fraction of the pulse period length. 
forwarding=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[Blanked_Train,Blank_Throughput] =Blank_it(Trains,Precapture_index,Pulse_Period,Blanking_window,Blanking_Shift,forwarding,Snips_fs);
figure
plot(Blanked_Train(:,:))
%% Clustering

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  n_clusters=2;
  emg_data= Trains;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STFT
        % Define STFT parameters.
        nperseg = 256;      % Length of each segment for STFT
        noverlap = 128;     % Overlap between segments
        stft_features_reshaped=STFT_Feature(emg_data,nperseg,noverlap,Snips_fs);  % 
        
%% K Mean clustering
        
        [cluster_labels,clustered_signals]=KCluster(stft_features_reshaped,emg_data,n_clusters);
        
%% DB indexing
        Cluster_ident=DB_Indexing(cluster_labels,emg_data,n_clusters);       

%% All trains in one plot stacked:
                if PDF_logic==1
                Color_pallet={'m', 'k' , 'b'};
                figure('units','inches','outerposition', [0 0 8.5 11])

                for i=1:1:size(Trains,1)
        
                        signal = Trains(i, :);
                        cluster_title=cluster_labels(i);
                        Cluster_color=Color_pallet{cluster_title};
                        stacked=normalize(signal,"range")+ i;

                        plot(stacked, 'DisplayName', sprintf('%s - #%d', cluster_title, i),'Color',Cluster_color);
                        hold on
                        text(length(stacked)+2,i,sprintf('# %g, Clsr %g',i,cluster_title))
                        hold on      
        
                end

                 title(sprintf('STACKED Ch%g Site%g Amp%g (part:%g) [%s]', Channel_Number,god,Amp_interest,Train_info(1,2),Descriptor));
                 xlabel('Index');
                 ylabel('Signal Amplitude');
                 grid on
                 add(Rep, Figure)
               
                 close all
                end

    


%% Data assignment for saving

        Mean_Cell{Amp_interest,Channel_Number}(god,:)= mean(clustered_signals{Cluster_ident},1);
        Median_Cell{Amp_interest,Channel_Number}(god,:)= median(clustered_signals{Cluster_ident},1);

        [RR,~]=find(cluster_labels==Cluster_ident);
        Selected_Train_number_cell{Amp_interest,Channel_Number}{god,1}=Train_info(RR,3); % this is absolute number

        Selected_Train_Matrix(Train_info(RR,3),Channel_Number)= 1;

         



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %end for bypass condition
                            end%end of site loops
        end % end of the Amp loop
                if PDF_logic==1
                close(Rep)  %closes the PDF file for the channel
                end
end %end of channel loop


%% Saving

Type='Cleared_Automatic';
        save(sprintf('AutoCleanedData/%s_AutoCleaned_(%s)_Raw%g.mat',Case,datestr(now, 'dd-mmm-yyyy'),Unblank_raw),'Mean_Cell','Median_Cell','Selected_Train_number_cell','Selected_Train_Matrix','Type','Descriptor','Snips_fs','Version_file','Unblank_raw')
        










%% STFT Clustering 
function stft_features_reshaped=STFT_Feature(emg_data,nperseg,noverlap,Snips_fs)

    stft_features = cell(size(emg_data, 1), 1); % Initialize an array to store STFT features for all samples.
    
    for i = 1:size(emg_data, 1) % Perform STFT analysis for all EMG signals and store the features.
        [S, F, T, P] = spectrogram(emg_data(i, :), hamming(nperseg), noverlap, nperseg, Snips_fs);
        stft_features{i} = abs(P);
    end
    % Convert the cell array of STFT features to a 3D matrix.
    num_samples = size(emg_data, 1);
    num_freq_bins = size(stft_features{1}, 1);
    num_time_bins = size(stft_features{1}, 2);
    stft_features_array = zeros(num_samples, num_freq_bins, num_time_bins);
    
    for i = 1:num_samples
        stft_features_array(i, :, :) = stft_features{i};
    end
    
    % Reshape the feature matrix for K-Means clustering.
    stft_features_reshaped = reshape(stft_features_array, num_samples, []);
end

%% K-means Clustering
function [cluster_labels,clustered_signals]=KCluster(stft_features_reshaped,emg_data,n_clusters)
    % Perform K-Means clustering on the STFT features.
    rng(0);  % Set the random seed for reproducibility
    opts = statset('Display', 'final', 'UseParallel', false);
    cluster_labels = kmeans(stft_features_reshaped, n_clusters,'Start', 'plus', 'Replicates', 5, 'Options', opts); 
    
    % Create a cell array to store cluster-specific signals.
    clustered_signals = cell(n_clusters, 1);
    
    % Populate the cell array with signals for each cluster.
    for cluster_id = 1:n_clusters
        cluster_signals = emg_data(cluster_labels == cluster_id, :);
        clustered_signals{cluster_id} = cluster_signals;
    end
end

%% DB indext identification

function   Cluster_ident_DB=DB_Indexing(cluster_labels,emg_data,n_clusters)
        %% Method to assess which cluster contains meaningful information after 
        % clustering is to use the "Davies-Bouldin Index" (DB Index). 
        % The DB Index measures the average similarity between each cluster and its most similar cluster, 
        % with a lower DB Index indicating better separation and, therefore, potentially more meaningful clusters. 
        %
        % Number of clusters (2 in this case).
        
        
        % Calculate the centroids for each cluster.
        cluster_centroids = zeros(n_clusters, size(emg_data, 2));
        
        for i = 1:n_clusters
            cluster_centroids(i, :) = mean(emg_data(cluster_labels == i, :));
        end
        
        % Initialize the DB Index to a high value.
        db_index = Inf;
        
        % Calculate the DB Index for each cluster.
        for i = 1:n_clusters
            for j = 1:n_clusters
                if i ~= j
                    % Calculate the average distance between the centroids of clusters i and j.
                    avg_distance = norm(cluster_centroids(i, :) - cluster_centroids(j, :));
        
                    % Calculate the average intra-cluster distance for cluster i.
                    % each pair of data point within a cluster vs the centroid of
                    % the cluster
                    avg_intra_distance_i = mean(pdist2(emg_data(cluster_labels == i, :), cluster_centroids(i, :)));
        
                    % Calculate the average intra-cluster distance for cluster j.
                    % each pair of data point within a cluster vs the centroid of
                    % the cluster
                    avg_intra_distance_j = mean(pdist2(emg_data(cluster_labels == j, :), cluster_centroids(j, :)));
        
                    % Calculate the Davies-Bouldin Index for cluster i.
                    db_i = (avg_intra_distance_i + avg_intra_distance_j) / avg_distance;
        
                    % Update the DB Index if cluster i has a lower DB Index than the previous clusters.
                    if db_i < db_index
                        db_index = db_i;
                        meaningful_cluster = i;
                        noise_cluster = j;
                    end
                end
            end
        end
        
           Cluster_ident_DB=meaningful_cluster;

end


% Blanking



function [Blanked_Train,Blank_Throughput] =Blank_it(Trains,Precapture_index,Pulse_Period,Blanking_window,Blanking_Shift,forwarding,Snips_fs)

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