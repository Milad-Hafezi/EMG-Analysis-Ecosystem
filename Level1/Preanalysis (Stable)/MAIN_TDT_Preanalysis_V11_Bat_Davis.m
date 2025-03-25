% there are some extra info that you can find in DATA/info
clear all

Num_parts=input('how many parts are your files? ');
initial_folder=pwd;

for J=1:1:Num_parts
    Directory{J}=uigetdir(pwd,sprintf('Please load part number %g',J));
    cd(Directory{J})
end

cd(initial_folder)

mkdir Raw
mkdir Rectified
mkdir Original

disp(Directory)

input('Press Enter to continue if the files are correct')
for Part=1:1:Num_parts
clearvars -except Num_parts Part Directory



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Channel_Numbers= 0             % Enter the specific channel that you want to look into. (FOR ALL CHANNELS TYPE 0)
% Unblank_raw= 0  % 0 for blanked and 1 for raw
% Part=0
% Part  =input('Part Number?   ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DATA = TDTbin2mat(Directory{Part}, 'TYPE', {'snips','epocs'}','CHANNEL',Channel_Numbers);
%DATA = TDTbin2mat(Directory, 'TYPE', {'snips','epocs'}');


     Case= DATA.info.Subject;
    %% Saving
    save(sprintf('Original/Matlab_(%s)_Part%g.mat',Case,Part),'-v7.3')
%file= uigetfile('*.mat');
%% Extracting th information Vectors
tic
for Unblank_raw=0:1:1
    load(sprintf('Original/Matlab_(%s)_Part%g.mat',Case,Part))
Count_Matrix (:,1)=DATA.epocs.CnA_.onset;
Count_Matrix (:,2)= DATA.epocs.CnA_.data;

Amp_Matrix(:,1)=DATA.epocs.AmpA.onset;
Amp_Matrix(:,2)= DATA.epocs.AmpA.data;


Note_codes= DATA.epocs.Note.data;
Note_data = DATA.epocs.Note.notes;
Note_time = DATA.epocs.Note.onset;

if Unblank_raw == 0
%Blanked Data
Snips_Times=DATA.snips.BMAT.ts;
Snips_fs = DATA.snips.BMAT.fs;
Snips_Channel_Code= DATA.snips.BMAT.chan;
Snips_Data = DATA.snips.BMAT.data;
elseif Unblank_raw==1
% Unblanked Data
Snips_Times=DATA.snips.uMAT.ts;
Snips_fs = DATA.snips.uMAT.fs;

Snips_Channel_Code= DATA.snips.uMAT.chan;
Snips_Data = DATA.snips.uMAT.data;
end

%Case number 
Case= DATA.info.Subject;


%% Breaking the big Snips matrix into 16 channel specific ones

for ch=1:1:16
    
CH=find(Snips_Channel_Code ==ch);
% UCH= find(UBSnips_Channel_Code ==ch);
Channel{ch}=Snips_Data(CH,:);
% UBChannel{ch}=UBSnips_Data(UCH,:);
end

%% Cleaning up the Snips time vector
Time=unique(Snips_Times);
% UBTime=unique(UBSnips_Times);
%% Epoc and Snip time difference compensation
Compensation=Time(1)-Amp_Matrix (1,1);
Time=Time-Compensation;

%% Detecting the site time range

Notes = string(Note_data);
Site_loc=find(Notes=='New Site');
Site_time= Note_time(Site_loc);
% Making sure that the difference between sites are correct, otherwise it
% means that the key was misspressed or been wrongly pressed. The
% difference should be at least one whole set of Stim precedure which is
% around 115 seconds
Difference= diff(Site_time);
Miss_pressed=find(Difference < 115)+1;
Site_time(Miss_pressed)=[];
% building the site range
Site_number=size(Site_time,1)
Site_Range=nan(Site_number,2);
Site_Range(:,1)=Site_time;
Site_Range(1:Site_number-1,2)=Site_time(2:Site_number,:);


%% Syncing the Data size
[Lia,Locb]=ismember(round(Amp_Matrix(:,1),2),round(Time,2)); %finding unique times based on the snips (it tell you that in Amp_Matrix, these X elemnts do exsist somewhere in the snip time)
[row,col]=find(Lia>0); % finding their location in Amp matrix
Amp=Amp_Matrix(row,2); %building the correct Amp Vector
Amp_Matrix_Revised=Amp_Matrix(row,:); % this builds an updated Matrix with time markers
Amp_Time=Amp_Matrix(row,1);
Count=Count_Matrix(row,2); %building the correct Count Vector
Count_Matrix_Revised=Count_Matrix(row,:);  % this builds an updated Matrix with time markers

% reverse Sync
[Liat,Locbt]=ismember(round(Time,2),round(Amp_Time(:,1),2)); % This tells you that which of the Time elements are not existant in the Amp matrix time
[rowt,colt]=find(Liat>0); % finding their location in Amp matrix
[rowtt,coltt]=find(Liat==0);
% Now I have to blank that train everywhere
Original_Time=Time;
Time(rowtt,:)=[];% Getting rid of the extra Times without any AMP/Count marker


for i=1:1:size(Channel,2) % Getting rid of the extra trains without any AMP/Count marker
    Channel{1, i}(rowtt,:)=[];

end

%% indexing the site positions for easier use in analysis

for sn=1:1:Site_number
   
    
    if sn<Site_number
        
    Loco=find(Time>Site_Range(sn,1)  & Time<Site_Range(sn,2));
    Site_index(sn,1)=min(Loco);
    Site_index(sn,2)=max(Loco);
    
    elseif sn==Site_number   
    
    Loco=find(Time>Site_Range(sn,1));
    if ~isempty(Loco)
    Site_index(sn,1)=min(Loco);
    Site_index(sn,2)=max(Loco);
    end
    end
    
end
%% finding the Snip_Range

   for bn=1:1:Site_number
       if bn == Site_number
           Ling = find(Time> Site_Range(bn,1));
       else
       Ling = find(Time> Site_Range(bn,1)  & Time < Site_Range(bn,2));
       end

       if ~isempty(Ling)
       Snip_Range(bn,1)=Time(Ling(1));
       Snip_Range(bn,2)=Time(Ling(end));
       Snip_Range(bn,3)=size(Ling,1);
       end
       
   end

%% Saving

clearvars -except file Unblank_raw rowtt Original_Time Case Channel UBChannel Time UBTime Site_Range Snip_Range General_Count General_Amp Site_index Site_number ...
    Count_Matrix Amp_Matrix Snips_fs UBSnips_fs Amp Count Compensation Amp_Matrix_Revised Count_Matrix_Revised Part Num_parts Directory
    if Unblank_raw==0  % saving unblanked
        save(sprintf('Rectified/Part%g SuperV10TDT (%s) %s.mat',Part,datestr(now, 'dd-mmm-yyyy'),Case),'-v7.3')
    elseif Unblank_raw ==1   %saving blanked
        save(sprintf('Raw/Part%g SuperV10TDT (%s) %s UNBLANKED_Raw.mat',Part,datestr(now, 'dd-mmm-yyyy'),Case),'-v7.3')
    end
end

end



toc

% Channel %all the trains are here
% Time % these are the times that each train is being delivered
% 
% Site_Range %Site time ranges
% Site_index %Site index ranges
% Site_number
% 
% Count %number of pulses in a train
% 
% Amp %Amp of each delivered train
% 
% Snips_fs % Recording frequency

