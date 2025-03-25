% there are some extra info that you can find in DATA/info
clear all
Address='D:\TDT\Synapse\Tanks\EMGDelayStimBatDevelopmentV6-230725-181642\Bat-test-2023-SFU-230728-011200';
B='Y'

prompt='Do you need Snips? Y/N  (Case sensitive):  ';
Check=input(prompt,'s');

if Check==B
    Channel_Numbers=input('Enter the specific channel that you want to look into. (FOR ALL CHANNELS TYPE 0):   ');
    DATA = TDTbin2mat(Address, 'TYPE', {'snips','epocs'}','CHANNEL',Channel_Numbers);
    %Case number 
    Case= DATA.info.Subject;
    %% Saving
    save(sprintf('%s.mat',Case),'-v7.3')
end

clear DATA
prompt='Do you need STREAMS? Y/N (Case sensitive):  ';
Check2=input(prompt);
T1=15000
T2=17000
if Check2==B

    data = TDTbin2mat(Address, 'STORE', 'fEMG', 'T1', T1, 'T2', T2);
     %Case number 
    Case= data.info.Subject;
    %% Saving
    save(sprintf('%s STREAM.mat',Case),'-v7.3')
end