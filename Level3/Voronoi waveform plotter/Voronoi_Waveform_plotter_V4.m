


tic

clearvars -except Meta_Data Descriptor Mean_Cell Median_Cell Selected_Train_Matrix Selected_Train_number_cell Snips_fs Type Version_file Unblank_raw

import mlreportgen.report.*
import mlreportgen.dom.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Butter_Cut_Off=45;   % remember that this filtration is only valid for this type of plots. This variable repeats throughout the code
Order=4;   
Scale_lining=1;
Normalized_Train=1;
Stacked_filtration=1;
    lowCutoff =45; % Lower cutoff frequency
    %highCutoff = 45;  % Upper cutoff frequency

Pre_filt=0;
    %lowCutoff_pre = 10; % Lower cutoff frequency
    highCutoff_pre = 20;  % Upper cutoff frequency

Attinuation=1;
Scale_map_ratio=[1 2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input('check Style on line 184 before you do anything else.')             %Style='TriMean';  %Do not forget to change this too


%% load the Meta Data file
if ~exist('Meta_Data')

[Partfile,Partpath] = uigetfile('*.mat','load the MetaData files');
                Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Channel_name);

end


% if ~exist('Version_file')
% [Partfile,Partpath] = uigetfile('*.mat','load the Cleaned up Data files');
%                 Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
%                 load(Load_Channel_name);
% 
% Version_file=Partfile;
% end


%% Find the amplitude levels and sites
%since it should be similar for all channels, I only considered the first
%channel

       Test  = Meta_Data{1, 1}.Guide_Matrix;
             Unique_Sites_MetaData= unique(Test(:,1));
             Unique_Amps_MetaData=unique(Test(:,6));
             Channels=size(Meta_Data,2);
             disp('Detected Amplitudes:  ')
             disp(Unique_Amps_MetaData)

%% Making directories
mkdir VoronoiWaveformMAPS
%% Run the waveform plotter
 Milad_Davis_code_UpgradeV7_Custom_VoronoiWaveform

%% Calculating the average equivalant square space to map the waveforms
clear Area
  for id = 1:size(Extract,1)
    Test_Area=polyshape([Extract{id, 6}(:,1)] ,[-Extract{id, 6}(:,2)]);
     clear area
     Area(id)=area(Test_Area);
  end

Equivalant_XY=sqrt(mode(Area));

%% Color pallet
    Grey_Back=  [0.909803921568627   0.905882352941176   0.898039215686275];
    Grey_lines=  [0.819607843137255   0.827450980392157   0.831372549019608];
    Grey_Unfair= [0.498039215686275   0.498039215686275   0.498039215686275];
    

                     
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
    
              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Here is the core code.
        for Amp_interest=256%Unique_Amps_MetaData' % remember that you need to have a horizontal vecotr to pull this off
                     [row,col, val] =find (Guide_Matrix(:,6)==Amp_interest);
                                      
                     clear Sub_matrix_Trains Sub_matrix_info
                     Sub_matrix_info= Guide_Matrix(row,:);
                     Sub_matrix_Trains=Guide_trains(row,:);

            %% plot the voronois
                 fig = figure;
                         hold on;
                         for id = 1:size(Extract,1)
                            %plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                            p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                            p.EdgeColor=Grey_lines;
                           %plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                           
                           %text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
                                                
                         end
                         if Scale_lining==1
                              line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
                              line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
                         end
                        axis off
                        axis equal
                            
                        xl = xlim; % getting the X and Y limits to put in for later plots
                        yl= ylim;
            
                        hold on
%% work from here on:

% You need to import the rest of the code and make is work with this
% It looks like the old code takes all the waveforms in One single matrix "Working_Train_set"
% So try to do the site loop first and then pic up wehre you left off only after the site loop is done and the appropreate measute (mean/meadian....) is done
% 

                            clear Trian_Site Train_info Verdic_site

                           for     god=Unique_Sites_MetaData'     % remember that you need to have a horizontal vecotr to pull this off    

                                                % Finding the approaved Trains
                                                clear Trains 

                                                     [rro,co,v]=find (Sub_matrix_info(:,1)==god);
                        
                                                    if isempty(rro)
                                                        disp('Amp not found for this site')
                                                    else
                                                        Train_info=Sub_matrix_info(rro,:);
                                                        Trains=Sub_matrix_Trains(rro,:);
                                                    end


%                                              


                                                if Pre_filt==1
                                                    clear Butter_Train
                                                      
                                                        % Design the bandstop filter
                                                        Wn=highCutoff_pre/(Snips_fs/2); % Filter parameters
                                                                [b,a] = butter(Order,Wn,'high');

                                                         for j=1:1:size(Trains,1)
                        
                                                                Butter_Train(j,:)= filtfilt(b,a,Trains(j,:));
                        
                                                          end
                                                            if exist('Butter_Train')
                                                                    Trains= Butter_Train;
                                                            end
    
                                                 end

                                                if isempty(Trains)
                                                Trian_Site(god,:)=zeros(1,size(Guide_trains,2));
                                                else

%                                                     coeff = pca(Trains,'Rows','pairwise');
%                                                     representativeSignal =trimmean(coeff',80,1); %mode (coeff',1); %trimmean(coeff',70,1); % Use the first principal component


                                                Trian_Site(god,:)= trimmean(Trains,40,1);  %geomean(Trains,1);   %representativeSignal;          % mode (Trains,1);          % trimmean(Trains,40,1);     %median(Trains,1);    
                                                end
                                                                           Style='TriMean';  %Do not forget to change this too
                        
                                                Verdic_site(god)=unique(Train_info(:,12));
                                                 
%                                               Verdict_Channel(Selected_Channel)=unique(Meta_Data{1, Selected_Channel}.Guide_Matrix(Selected_Trians{1,1},13));
                          
                            end%end of site loops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%% Filtration for all selected Train_sites
                                  if Stacked_filtration==1
                                                   clear  Butter_Train
%                                                        d = designfilt("bandpassfir", ...
%                                                                         FilterOrder=Order,CutoffFrequency1=lowCutoff, ...
%                                                                         CutoffFrequency2=highCutoff,SampleRate=Snips_fs);
                                                                Wn=lowCutoff/(Snips_fs/2); % Filter parameters
                                                                [b,a] = butter(Order,Wn); % Set as butterworth filter
                                                            for j=1:1:size(Trian_Site,1)
                        
                                                                Butter_Train(j,:)= filtfilt(b,a,Trian_Site(j,:));
                        
                                                            end
                                                            Trian_Site = Butter_Train;
                                 end
%% plot the waveforms


            clear Calculated_time  Working_Train_set Normalized_Calculated_time
        
            [Working_Train_set, Cen,Scale]=normalize(Trian_Site,2,"range",[0 round(Equivalant_XY)]); %When A is an array, normalize returns C and S as arrays such that N = (A - C) ./ S. Each value in C is the centering value used to perform the normalization along the specified dimension. For example, if A is a 10-by-10 matrix of data and normalize operates along the first dimension, then C is a 1-by-10 vector containing the centering value for each column in A.
            %... [N,C,S] = normalize(___)  where   N = (A - C) ./ S
            %   OR    n.*s+c
            
            Relative_Scale=normalize(Scale,'range',Scale_map_ratio);
            

                                  
            Calculated_time=1:1:size(Working_Train_set,2);
            Normalized_Calculated_time=normalize(Calculated_time,2,"range",1.4*[0 round(Equivalant_XY)]);
            % Calculating the on and off index for the stim box
                    size(Working_Train_set,2)
                    start_index=round(0.5* (Snips_fs)); % start index for stim
                    end_index=round(1*(Snips_fs));  % end of stim index
    for J=1:1:2
              for i=1:1:size(Working_Train_set,1)
                    clear ro
                    [ro,val]=find(cell2mat(Extract(:,7))==i);
                        if isempty(ro)
                            
                            for Ka=1:1:size(Extract(:,8),1)
                                if isempty(ro)
                                [ro,val]=find(cell2mat(Extract(Ka,8))==i);
                                    if ~isempty(ro)
                                        found=Ka;
                                    end
                                end
                            end

                            ro=found;
                        end

                        if isempty(ro)
                            disp('PROBLEM DETECTED S0012')
                        end


                        if Attinuation==0
                            Scaler=Relative_Scale(i);
                        else
                            if Verdic_site(i)==0
                                Scaler=Scale_map_ratio(1)/10;
                            else
                                Scaler=Relative_Scale(i);
                            end
                        end
                     X_interest=Extract{ro,2};
                     Y_interest=-Extract{ro,3};
                     Shifted_X=Normalized_Calculated_time+X_interest-Normalized_Calculated_time(start_index);
                     Shifted_Y=(Scaler*Working_Train_set(i,:))+Y_interest;
                        % stim box

                      Rec_position=[X_interest, min(Shifted_Y), Normalized_Calculated_time(end_index)- Normalized_Calculated_time(start_index),  max(Shifted_Y)-min(Shifted_Y)]; % this is the height



%                        rectangle('Position',Rec_position,'FaceColor', ([227 229 229]/256),'EdgeColor','none')
                   if J==1
                   rectangle('Position',Rec_position,'FaceColor', ([220 220 220]/256) ,'EdgeColor','none')
                  
                   hold on
                   elseif J==2
                        if Stacked_filtration==1
                            Line_thickness=0.2;
                        else
                            Line_thickness= 0.01;
                        end

                    plot(Shifted_X,Shifted_Y, "LineWidth",Line_thickness)
                    hold on
                     

                   end
              end
    end
 
    %% 
    for id = 1:size(Extract,1)
       
%     plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
        
       
        %% check for positive redo sites
        

        if ~Extract{id,9}==0
            
            clear Redo_thresh  
            for i=1:1:size(Extract{id,9},2)

                [rR,~]=find(Guide_Matrix(:,1)==Extract{id,9}(i));
                      
                Redo_thresh(i)=Guide_Matrix(rR(1),13); % C{Selected_Channel}(Extract{id,9}(i),4);  
            end

                [roo,cool,v] = find(Redo_thresh>0 & Redo_thresh<=Amp_interest);


                        if nnz(roo) > 0
                             [rR,~]=find(Guide_Matrix(:,1)==Extract{id,7});

                            Activity= Guide_Matrix(rR(1),12); %C{Selected_Channel}(Extract{id,7},3);
                        else
                            Activity=0;
                        end
                        
        else % if there is no redo we just find weather the site is active or if it has the correct threshold
                     [rR,~]=find(Guide_Matrix(:,1)==Extract{id,7});

                   Verdict_threhs=Guide_Matrix(rR(1),13); %C{Selected_Channel}(Extract{id,7},4);
                   Act_check= Amp_interest-Verdict_threhs;
                   if Act_check < 0
                       Activity=0; %negative site for my amplitude
                   elseif Act_check >= 0
                       Activity= Guide_Matrix(rR(1),12);    %C{Selected_Channel}(Extract{id,7},3);
                   end
        end
       


      if   Activity ==1
          Activity_marker='*-';
      else
          Activity_marker='';
      end



      text((Extract{id,2}+1),-(Extract{id,3}+(0.10*Equivalant_XY)),sprintf('%s %s',Activity_marker,Extract{id,1}),"FontSize",1,Color=Grey_Unfair) 

    end





    %% save gfc
        
            % making sure the folder is created
            fullPath = fullfile('VoronoiWaveformMAPS', Style);

            mkdir(fullPath)
  
            title(sprintf('Raw%g %s Waveforms (%s) uA:%g Ch%g Filter=%g Prefilt=%g \n  %s',Unblank_raw,Style,Case,Amp_interest, Channel_Number,Stacked_filtration,Pre_filt,Muscle_name),'FontSize',5)
            set(gcf, 'Renderer', 'painters');  % Vector graphics renderer (better for vector formats)
            
            saveas(gca,sprintf('VoronoiWaveformMAPS/%s/Raw%g_(%s)Woronoiwaveform uA%g Ch%g (%s)Filter%g Pre%g.svg',Unblank_raw,Style,Case,Amp_interest, Channel_Number,Muscle_name,Stacked_filtration,Pre_filt)) %  SVG
       


close all


        end % end of the Amp loop
end %end of channel loop


Style

toc













