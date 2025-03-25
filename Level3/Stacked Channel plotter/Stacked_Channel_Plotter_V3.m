
%  This code selectively plots the all channels per selected site.


tic

clearvars -except Meta_Data Descriptor Mean_Cell Median_Cell Selected_Train_Matrix Selected_Train_number_cell Snips_fs Type Version_file Unblank_raw
close all
import mlreportgen.report.*
import mlreportgen.dom.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Separation_type='Tec'; % seperates trains based on the max height of the previouse train
Separation=1;          % Ratio of the seperation  (this multiplies with the max height determined by the previouse train)
Stacked_filtration=0; % prefiltration
Post_filt=0;   %post filtration
Butter_Cut_Off=45;   % remember that this filtration is only valid for this type of plots. This variable repeats throughout the code
Order=4;   

Normalized_Train=1;
PDF_logic=0; % Tunrs the PDF report generation On (1) or off (anything but 1)
% More info on PDF Reports
% https://www.mathworks.com/help/rptgen/ug/center_figure_snapshot_on_a_page.html
% https://www.mathworks.com/help/rptgen/ug/mlreportgen.dom.pdfpagelayout-class.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% load the Meta Data file
if ~exist('Meta_Data')

[Partfile,Partpath] = uigetfile('*.mat','load the MetaData files');
                Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Channel_name);

end

% 
% if ~exist('Version_file')
% [Partfile,Partpath] = uigetfile('*.mat','load the Cleaned up Data files');
%                 Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
%                 load(Load_Channel_name);
% Version_file=Partfile;
% end

if ~exist('Version_file')
    Type='U';
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

%% Making directories
mkdir Stacked_muscles
mkdir Stacked_musclesSVG

%%   User inputs
                    prompt='Which Site are you plotting the stacked plots for? ';
                    god=input(prompt);
                    disp(Unique_Amps)
                    Amp_interest=input('What amplitude Should I use?  ');
                   
%% PDF Report setup
                if PDF_logic==1
                Rep = Report(sprintf('Stacked_muscles/StackedCH_%s_%s_Amp_%g_Site%g_Norm%g_Postfilt%g',Case,Type,Amp_interest,god ,Normalized_Train,Post_filt),'pdf');  % add the apmplitude


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
                 %% finding the Amplitude for that channel


                 [row,col, val] =find (Guide_Matrix(:,6)==Amp_interest);
                                      
                     clear Sub_matrix_Trains Sub_matrix_info
                     Sub_matrix_info= Guide_Matrix(row,:);
                     Sub_matrix_Trains=Guide_trains(row,:);
                        
            
                     %% finding the correct site
                            
                     clear Trains   Train_info
                                            % finding the correct sites
                                               [ro,co,v]=find (Sub_matrix_info(:,1)==god);
                        
                                                    if isempty(ro)
                                                        disp('Amp not found for this site')
                                                    else
                                                        Train_info=Sub_matrix_info(ro,:);
                                                        Trains=Sub_matrix_Trains(ro,:);
                                                    end





                        %% Finding the approaved Trains (Using Cleaned Data
%                          clear Trains 
%                          Selected_Trians=Selected_Train_number_cell{Amp_interest,Selected_Channel}(god,1);
%                          Trains=Meta_Data{1, Selected_Channel}.Guide_trains(Selected_Trians{1,1},:);


%%
                         if Stacked_filtration==1
                                 Wn=Butter_Cut_Off/(Snips_fs/2); % Filter parameters
                                [b,a] = butter(Order,Wn); % Set as butterworth filter

                                    for j=1:1:size(Trains,1)

                                        Butter_Train(j,:)= filtfilt(b,a,Trains(j,:));

                                    end
                                    Trains = Butter_Train;
                         end

                        Trian_Channel(Selected_Channel,:)= trimmean(Trains,40,1);  %median(Trains,1);
                                                   Style='Trimmean';  %Do not forget to change this too


                        Verdict_Channel(Selected_Channel)=unique(Train_info(:,13));
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
end %end of channel loop
%% post filtration


        if Post_filt==1
            clear Butter_Train
             Wn=Butter_Cut_Off/(Snips_fs/2); % Filter parameters
             [b,a] = butter(Order,Wn); % Set as butterworth filter
                
              for j=1:1:size(Trian_Channel,1)

                                        Butter_Train(j,:)= filtfilt(b,a,Trian_Channel(j,:));

               end
            Trian_Channel=Butter_Train;
        end

%% plotting in one shot

 %% variable assignment
              
                 Site_No=god;
                   if Normalized_Train==1
                        Train= normalize(Trian_Channel,2,"range");
                   end

                   
                 %%  Stacked plots
                
                    if strcmp(Separation_type , 'Tec')
                        Height_vector=abs(max(Train,[],2)-min(Train,[],2));
                        Vertical_Separation_vector= Separation*Height_vector;
                  
                    end
                
                    clear Separated_Train  Seperated_text Time_train
                    Time_train=-500:(1/Snips_fs)*1000:1000; % time train in ms
                    Separated_Train=zeros(size(Train));
                    figure('units','inches','outerposition', [0 0 8.5 11])
                  
                
                    for no= 1:1:size(Train,1) 
                        if no==1
                            Separated_Train(1,:)= Train(1,:);
                            Seperated_text(1)= median(Train(no,:));
                
                        else
                            Separated_Train(no,:)= Train(no,:)+ sum(Vertical_Separation_vector(1:no-1));
                            Seperated_text(no)= median(Train(no,:))+ sum(Vertical_Separation_vector(1:no-1));
                        end
                    end
                    
                    plot(Time_train,Separated_Train,'LineWidth',0.5);
           
                    hold on
                    xline(0,'--r');
                    xline(500,'--r');
                    
                    title(sprintf('%s_%s_(%s) Stacked channels (Site%g), uA:%g,n=%g, Normalized:%g,PostFilt=%g-Raw%g', Style,Type,Case,Site_No,Amp_interest,size(Train,1),Normalized_Train,Post_filt,Unblank_raw));
                    xlabel('Time Relative to Stim onset (ms)');
                    ylabel('EMG (a.u)');
                    grid on
                    
                    for sop=1:1:size(Train,1)
                       
                      
                            text(1001,(Seperated_text(sop)),sprintf('(Ch%g) - %g',sop,Verdict_Channel(sop)),'Color','r','FontSize', 6)
                           
                    end



% Save SVG
                  saveas(gca,sprintf('Stacked_musclesSVG/%s_%s_(%s)_Stacked channels_(Site%g)uA%g_Norm%g,Prefilt%g_PostFilt%g_Raw%g.svg', Style,Type,Case,Site_No,Amp_interest,Normalized_Train,Stacked_filtration,Post_filt,Unblank_raw)) %exporting to svg



if PDF_logic==1
%report
add(Rep, Figure)

close(Rep)%closes the PDF file for the channel
end
close all


toc


