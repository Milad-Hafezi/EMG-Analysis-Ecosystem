
% 

tic

clear all

%% load the Verdict    (I have disabled this becasue of excel problems) make sure to put it back

[file,path] = uigetfile('*.xlsx','Select the Verdict file for your case');
cd(  path);
Verdict=importdata(file);

MUSCLES=Verdict.textdata.Sheet17(1:16,1);
MUSCLES_Abb=Verdict.textdata.Sheet17(1:16,3);
% Case=Verdict.textdata.Sheet17{1, 2} 
% Color=Verdict.data.Sheet17(1:16, 2:4) 
% Color=Color./256;

Check = struct2cell(Verdict); % just to know if there is any text with the verdict
if size(Check,1)>1    %seperating the text from numerical data
Verdict_Text=Verdict.textdata;
Verdict=Verdict.data;
end
C = struct2cell(Verdict); % transforming the Struct to cell
C_text=struct2cell(Verdict_Text);


if size(C,1)>16 % if Verdict has 17 sheets we only need 16 channels of it
   C(17)=[];
end



for Bin=1:1:size(C,1)
    clear Sanity
    Sanity=C{Bin}(:,4)<888;
    Coactivation(:,Bin)=C{Bin}(:,3).*Sanity;
    CoThresh(:,Bin)=C{Bin}(:,4);
end
%% clean up
clearvars -except Unblank_raw Original_Time rowtt Amp Case Channel Count Site_index Site_number Site_Range Snips_fs Time UBChannel UBSnips_fs UBTime Compensation Count_Matrix Amp_Matrix Amp_Matrix_Revised Count_Matrix_Revised ...
    Case_no Train_count Section_no CoThresh Coactivation MUSCLES_Abb  MUSCLES Guide_Matrix Guide_trains Guide_trains_downsampled  Meta_Data Unblank_raw
%%

mkdir MetaFiles
mkdir MetaFiles_Combined



%% Super loop for different parts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    prompt='How many file parts are there?   ';
    Part_numbers=input(prompt); %  first bat= 7
  
    Repeating=input('Do you want to star with fresh parts (enter 0) or load an exsisting directory address (press any keys)       ' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % inital working folder:
        Current_folder=pwd;


    if   Repeating==0
              clear Part_address_name
            for parts=1:1:Part_numbers
                [Partfile,Partpath] = uigetfile('*.mat',sprintf('Load Data file for part %g', parts));
                cd(Partpath);
                Part_address_name{parts}=sprintf('%s\%s',Partpath,Partfile);
                stat= strfind(Partfile,'Raw');
            end
            % saving the file data
            cd(Current_folder)
            if isempty(stat)
                raw=0;
            else
                raw=1;
            end
            save(sprintf("Part_address_name_Raw%g.mat",raw),"Part_address_name")
    else
        [Partfile,Partpath] = uigetfile('*.mat','load: Part_address_name');
                Load_Part_address_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Part_address_name);
                
    end

cd(Current_folder)

%% Batch vs Single
Batch= 0; % set 0 for single or anything else for batch running

%% inputs
% Constants & Inputs 1
if Batch==0
    prompt='Which channel to START from?';
    Channel_Start=input(prompt);
    prompt='Which channel to END with?';
    Channel_End=input(prompt);
    prompt='START site?';
    Start_Site=input(prompt);
    prompt='Enter END site or just leave it for full run: ';
    checks=input(prompt);
   

    Amp_Start=input('Starting Amplitude (Leave it for a full run)?  ');
    Amp_end= input('End amplitude(Leave it for a full run)?   ');
else
    prompt='Which channel do you want to analyze?';
    Channel_Start=input(prompt);
    Channel_End=Channel_Start;
    prompt='Which Site should I START from?';
    Start_Site=input(prompt);
    prompt='Which Site should I END with?';
    End_Site=input(prompt);
end
%% Initial allocaitions decide where to put these:

Train_count=0;
Full_notch =0;
Part_last_site=0;




for GOD=Channel_Start:1:Channel_End
%% Clearing the inside variables
 clearvars -except Unblank_raw Amp Case Channel Count Site_index Site_number Site_Range Snips_fs Time UBChannel UBTime UBSnips_fs ...
     GOD Channel_Start Channel_End Start_Site End_Site Rep ...
     Compensation Count_Matrix Amp_Matrix Amp_Matrix_Revised Count_Matrix_Revised...
     Amp_Start Amp_end Full_notch Case_no Section_no CoThresh Coactivation  MUSCLES MUSCLES_Abb Guide_Matrix Guide_trains Guide_trains_downsampled ...
     Part_last_site Partition Part_address_name checks Meta_Data Unblank_raw
     %Checkers CH interval  EMG_number  LOG Sync AMP Case_Number EMG_Channel_Numbers GOD Start_Site
 close all
 
Train_count=0;
Part_last_site=0;


 %% Super loop in parts STARTS HERE
for Partition=1:1:size(Part_address_name,2)  

%% Load the files
load(Part_address_name{Partition})

 if nnz(checks) == 0
        End_Site=Site_number;   
 else
        End_Site=checks;
 end
%Case_no=str2num(sprintf('%s%s',Case(1:2),Case(4:5)));
%%  Data assignement

EMG=double(Channel{GOD});
Channel_No= sprintf('Case_ %s_Ch%g ',Case,GOD); % the site number


%% Report Generation
% More info
% https://www.mathworks.com/help/rptgen/ug/center_figure_snapshot_on_a_page.html
% https://www.mathworks.com/help/rptgen/ug/mlreportgen.dom.pdfpagelayout-class.html
        if Full_notch ==1
            Notch_stat='NOTCHED';
        else
            Notch_stat='';
        end


%%  Inner loop (god)

for god=Start_Site:1:End_Site % Inner loop (site based)

clearvars -except Unblank_raw Amp Case Channel Count Site_index Site_number Site_Range Snips_fs Time UBChannel UBTime UBSnips_fs...
     GOD Channel_Start Channel_End Start_Site End_Site Rep ...
     EMG Channel_No...
     god...
     Compensation Count_Matrix Amp_Matrix Amp_Matrix_Revised Count_Matrix_Revised...
     Amp_Start Amp_end Full_notch Case_no Train_count Section_no CoThresh Coactivation MUSCLES MUSCLES_Abb Guide_Matrix Guide_trains Guide_trains_downsampled ...
     Part_last_site Partition Part_address_name checks Meta_Data Unblank_raw


     
% Defining Site Number   
Site_No=sprintf('Site%g',god);

    %% Constants & Inputs %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Developer_Mode='N';   % Either Y for yes or N for no
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spont_Detection =0;   % 1 enables the detection and elimination of train with spontaneous activity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Filter=1;             %Turns mean and STD filtaration on or off (1 for filteration after averaging and 2 for filteration before averaging)
Butter_Cut_Off = 45; %ORIGINALLY 180!   %Butterworth filteration cutt of frequency (low pass)
Order=4; % Butterworth filter order

Stacked_butter=0;% this would make it possible to filter each train in the stacked plot.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 60Hz cuttoffs
High_cut=123; % high pass
Low_Cut= 117;% low pass
Order_60=2 ; % order of the butterworth filter for filtering 60 hz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Precapture = 500; % Precapture of TDT in miliseconds (ms)
duration= 1500;  % Duration of Capture of TDT in miliseconds (ms)
Stim_window_time= 500; %  Stim window  (ms)
Section_times=100; % The length of the bins for sectioning (ms)

Window1=5;   % duration of the first window of analysis (ms)
Window2_start=20; % start of the 2nd window of analysis (ms)



Separation=25;%was 10 originally   %Separatation ratio based on the Max or STD value of EMG
Separation_STIM=4; %  Separatation ratio based on the STD during stimulation 
Mean_seperation_Ratio=2; % Seperation betwen Mean plots at the end

Separation_type= 'CON';   % if Maximume based: write 'Max' otherwise write 'STD'
                          % in the case of 'Max' do not forget to change
                          % the separation amount! 'STD' is from the backgro
                          % bacground noise. 
                          % 'CON' Seperates the trains based on the Std of
                          % the stim period (concatenated across different
                          % trials)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ylimit=1; % fixes the y limit on plots for the max amount of all channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paper=0; % if you have it at 1, it will plot the paper style figure and anything else doesn't 
Stacked_Plot =1; % plots the stacked plots for 1 and not for anything else.
Notching=0; % 1 to notch filter only the stacked mean plot and 0 to not to filter at all.
Stacked_mean_Plots=1;  % 1 to turn on any Stacked mean plot (the plot at the end) or 0 to turn it off
Downsample_rate=0;  % put 0 for no downsampling and any positive value for turning the downsampling on
Meta_Downsample_rate=5; % this controls the downsampling rate for the Meta data genaration.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Blanking
Blanking_control=0;
Pulse_Period=4.99712;   %ms     0.4897  0.00499712
Blanking_window=3.1; % ms
Blanking_Shift=43; % a number from 0 to 100 (percentage) this can shift all the blanking windows by a fraction of the pulse period length. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Assignment
% time-based transformation


Index_Start=Site_index(god,1);
Index_End=Site_index(god,2);
EMG_Site = EMG(Index_Start:Index_End,:); % this is the EMG trains for a specific site
Time_Tick= Time(Index_Start:Index_End);
Amp_Site= Amp(Index_Start:Index_End); %  Amplitudes of each train

Time_Start=Site_Range(god,1); % finding the correct amps by searching for their times
Time_End=Site_Range(god,2);

        if god == Site_number
           [Ling,lung] = find(Amp_Matrix_Revised(:,1) > Site_Range(god,1));
        else
           [Ling, lung] = find(Amp_Matrix_Revised(:,1)> Site_Range(god,1)  & Amp_Matrix_Revised(:,1) < Site_Range(god,2));
       end
       

Amp_Site_Matrix=Amp_Matrix_Revised(Ling,:); % this is the best way to look for a specific train or specific amplitude


%Amp_Matrix_Site=Amp_Matrix_Revised(Index_Start:Index_End)
Count_Site =Count(Index_Start:Index_End); % Number of pulses of that train
interval=1000/(Snips_fs);
Time_Train=-Precapture:interval:(duration-Precapture);

NAME='';

%% Detrending and Rectification
%EMG_Site=abs(detrend(EMG_Site));
 %EMG_Site=abs(EMG_Site);
%% Detecting the same amps and their max
    %UNI=flip(unique(Amp_Site))
    UNI=flip(sort(unique(Amp_Site_Matrix(:,2))));
    Max_UNI=max(UNI);
    
    %% finding the max of the signal
    Max=max(max(EMG_Site));
    
    %% the loop for each unique value stars here.
    %% finding the amplitude range

    if isempty(Amp_Start) || isempty(Amp_end)
        cil=1;
        cil2=length(UNI);
    else
    
     [cil man]=find(UNI==Amp_Start);
      [cil2 man2]= find(UNI==Amp_end);
    end

    for uni=cil:1:cil2%1:1:length(UNI)


        Site_Amp=sprintf('%g uA',UNI(uni)); %for plots and reports
    %% Finding the positions of the each Amp
         clear Train AM Train_pre Train_aft RMS_TRAIN RMS_TRAN_Pre 
         AM=find (Amp_Site_Matrix(:,2)==UNI(uni));
            Timies=Amp_Site_Matrix(AM,1); % Taking the Epoc times from Amplitude matrix
            for vi=1:1:length(Timies)
                clear BM
                BM=find(round(Time_Tick,2)==round(Timies(vi),2)); %finds the sister times in the Time vector
%                     if length(BM) > 1
%                         SANITY (vi)= 'NO';
%                     elseif length(BM) == 1
%                         SANITY(vi) = 'YES';
%                     else
%                         SANITY(vi) = 'NO';
%                     end
                if length(BM)>= 1
                 TICK(vi,:)=Time_Tick(BM);
                 Train(vi,:)=EMG_Site(BM,:);
                end

            end


    %% EMG Time Periods
        Precapture_index= round((Precapture/1000)*Snips_fs);
        Window1_Index=round((Window1/1000)*Snips_fs);
        Window2_start_Index=round((Window2_start/1000)*Snips_fs);
        duration_index= round((duration/1000)*Snips_fs);
        Stim_window_index = round((Stim_window_time/1000)*Snips_fs);

        Sectioning_Index=floor((Section_times/1000)*Snips_fs); % Using this to make a 100 X 100 Matrix


         
    %% Train Mass filtraion
        % Filtering out the possible 60 Hz and the harmonics for our data
        clear Notched_Train
        

        % detemining the strategy
        
        if Notching==1 || Full_notch ==1

                freq = Snips_fs; %  sampling frequency (Hz)
        
            
        Fnotch = 120;  % Notch Frequency
        BW     = 12;   % Bandwidth
        Apass  = 1.2;    % Bandwidth Attenuation

                  % [b, a] = iirnotch(Fnotch/(freq/2), BW/(freq/2), Apass);
                    
                 
                     W1=Low_Cut/(freq/2);
                     W2=High_cut/(freq/2);
                     WN=[W1 W2];
                     [b,a] = butter(Order_60,WN,'stop');
                    Notched_Train=zeros(size(Train));
                for filty=1:1:size(Train,1)
                    Notched_Train(filty,:)= filtfilt(b,a,Train(filty,:));
                end
% 
        end

    
        if Full_notch ==1
            Notching=0;
            Train=Notched_Train;
            Notch_stat='NOTCHED';
        else
            Notch_stat='';
        end

%% Blanking Trains 
% We try to blank the trains just in case there residual artifact. I do
% this by blanking a preset pattern starting from the stimulation onset.
if Blanking_control==1
% Logic vector:

% Blank_Logic= zeros(1,size(Train,2));
% finding zero index
Blanking_onset_index=(Precapture/1000)*Snips_fs; %ceil
Blank_Steps=(Blanking_window/1000)*Snips_fs; %floor
Blank_Throughput=ceil(((Pulse_Period-Blanking_window)/1000)*Snips_fs);
Window=(Pulse_Period/1000)*Snips_fs;
forwarding=0;
Blanking_Shift=(Blanking_Shift/100)*Window; %redefining the blanking_shift in index domain
Blanked_Train=Train;

for bil=1:1:100 % make sure to change this using Count matrix later if you are doing a short train stim
        Starting=floor(Blanking_Shift +Blanking_onset_index+forwarding); 
        Ending= ceil(Blanking_Shift+Blanking_onset_index+forwarding+Blank_Steps);
            for rows=1:1:size(Train,1)
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
Train=Blanked_Train;

end

% Deep filtration (this would filter each train)
clear Butter_Train
if Stacked_butter ==1
    Butter_Train=zeros(size(Train));
    freq = Snips_fs; %  sampling frequency (Hz)
    Wn=Butter_Cut_Off/(freq/2); % Filter parameters
    [b,a] = butter(Order,Wn); % Set as butterworth filter

    for jiz=1:1:size(Train,1)
            Butter_Train(jiz,:)= filtfilt(b,a,Train(jiz,:));
    end
            Train=Butter_Train;
end       

    %% Precapture baseline Sectioning, Mean and STD
        Train_pre=Train(:,1:Precapture_index);
    
            Baseline=reshape(Train_pre, [], 1);       
            Baseline_Mean= mean(Baseline);
            Baseline_STD= std(Baseline);
              % Prepration for ploting
            at=[Time_Train(1) Time_Train(end)];
            a_stdhigh =[Baseline_Mean+Baseline_STD,Baseline_Mean+Baseline_STD];
            a_stdlow = [Baseline_Mean-Baseline_STD,Baseline_Mean-Baseline_STD];
            a_mean= [Baseline_Mean Baseline_Mean];
     %% Stimulation period sectioning, mean and STD
          Train_During=Train(:,(Precapture_index+1):(Precapture_index+Stim_window_index));
    
            During_Stim=reshape(Train_During, [], 1);       
            Mean_During_Stim= mean(During_Stim);
            STD_During_Stim= std(During_Stim);
              % Prepration for ploting
            att=[Time_Train(1) Time_Train(end)];
            a_stdhigh_During =[Mean_During_Stim+STD_During_Stim,Mean_During_Stim+STD_During_Stim];
            a_stdlow_During = [Mean_During_Stim-STD_During_Stim,Mean_During_Stim-STD_During_Stim];
            a_mean_During= [Mean_During_Stim Mean_During_Stim];

     %% Sectioning #2 The stim and post
        Train_aft=Train(:,Precapture_index+Window2_start_Index:Precapture_index+Stim_window_index); % from the onset of stim to the end of it
    
 %% Averaging
clear Mean STD  Notched_Mean
 if size(Train,1) > 1 %%CHANGE HERE
     Mean= mean(Train); %Mean of each row (between trains)
     STD=std(Train,0,1); %STD of each column (Between Trains)
        if Notching==1
            Notched_Mean= mean(Notched_Train);
        end

 elseif size(Train,1)== 1
     
     Mean= Train;
     STD=std(Train,0,1);
        if Notching==1
            Notched_Mean= Notched_Train;
        end
 end
 %
 if Developer_Mode =='Y' %Developer Mode
    figure
    PLOT=plot(Time_Train,Mean,Time_Train,Train);
     PLOT(1).LineWidth = 2;
    title(sprintf('(%s)Superimposed trains vs Mean', NAME));
    xlabel('Time Relative to Stim onset (ms)');
    ylabel('EMG & Mean (uV)');
    
    grid on

    figure
    plot(Time_Train,Mean);
    title(sprintf('(%s)Mean',NAME));
    xlabel('Time Relative to Stim onset (ms)');
    ylabel('Mean (uV)');
    grid on
 end
 
 %% Mean Filtration and butterworth filter setup
freq = Snips_fs; %  sampling frequency (Hz)
Wn=Butter_Cut_Off/(freq/2); % Filter parameters
if Filter==1
[b,a] = butter(Order,Wn); % Set as butterworth filter

    Filtered_Mean= filtfilt(b,a,Mean);
    Filtered_STD= filtfilt(b,a,STD);

        if Notching==1
            Notched_Mean_filtered= filtfilt(b,a,Notched_Mean);
        end

% 
  if Developer_Mode =='Y' %Developer Mode
    figure
    plot(Time_Train,Filtered_Mean);
    title(sprintf('(%s)Filtered Mean', NAME));
    xlabel('Time Relative to Stim onset (ms)');
    ylabel('Filtered Mean EMG (uV)');
    grid on
  end
else
  Filtered_Mean=Mean;
  
end




%% logging Meta Data

clear Section_no
Section_no=GOD*(10^6)+god*(10^3)+UNI(uni)
for Chopp=1:1:size(Train,1)
    clear Vector
    Train_count=Train_count+1;

        if CoThresh(Part_last_site+god,GOD) > UNI(uni) % making sure that the verdict is  sperated for  the subthreshold trains vs the superthrehshold and thrshold trains
            VerdictYN=0;
        else
            VerdictYN=Coactivation(Part_last_site+god,GOD);
        end

    Vector=[Part_last_site+god, Partition,	Train_count,	GOD, 	god, 	UNI(uni),	TICK(Chopp),	Chopp,	Precapture_index,	(Precapture_index+Stim_window_index),	Section_no, 	VerdictYN, 	CoThresh(Part_last_site+god,GOD), Meta_Downsample_rate, 0,0,0];  % please note that the zeros are place holders

    Guide_Matrix(Train_count,:)=Vector;

    Guide_trains(Train_count,:)=Train(Chopp,:);
%     Guide_trains_downsampled(Train_count,:)=  downsample(Train(Chopp,:),Meta_Downsample_rate) ;
   

end



end % end of the uni (Amplitude) loop





end % end of god (site) loop
Part_last_site=Part_last_site+god; % compiling the number of real sites.
end  %super loop ends here  FOR EACH PART

    Muscle_name=MUSCLES{GOD};
    Channel_Number=GOD;
    save(sprintf('MetaFiles/Meta_info_Matrix_CH(%g)_(%s)_%s_Raw%g.mat',GOD,datestr(now, 'dd-mmm-yyyy'),Case,Unblank_raw),'Guide_trains','Guide_Matrix','Muscle_name','Case','Channel_Number','Snips_fs','Unblank_raw','-v7.3')
%    save(sprintf('MetaFiles/Meta_Train_CH(%g)_(%s)_%s.mat',GOD,datestr(now, 'dd-mmm-yyyy'),Case),'Guide_trains','Muscle_name','Case','-v7.3')
%     save(sprintf('MetaFiles/Meta_Downsampled_Train_CH(%g)_(%s)_%s.mat',GOD,datestr(now, 'dd-mmm-yyyy'),Case),'Guide_trains_downsampled','Muscle_name','Case','-v7.3')


    
    
Struct_Vector=struct('Guide_trains',Guide_trains,'Guide_Matrix',Guide_Matrix,'Muscle_name',Muscle_name,'Case',Case, 'Channel_Number',Channel_Number,'Snips_fs',Snips_fs);
Meta_Data{Channel_Number}=Struct_Vector;

 



    clear Guide_Matrix Guide_trains...
           Struct_Vector   %Guide_trains_downsampled
    clc





end % end of GOD (channel) loop
    save(sprintf('MetaFiles_Combined/Meta_Data(%s)_%s_Raw%g.mat',Case,datestr(now, 'dd-mmm-yyyy'),Unblank_raw),'Meta_Data','Unblank_raw','-v7.3')

toc
%% TO DO

% save each meta file seperately for a channel
% Have an absolute site number


