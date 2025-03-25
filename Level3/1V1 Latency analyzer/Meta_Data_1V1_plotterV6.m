%% Initialize
clearvars -except Meta_Data
close all

%% load the meta data
if exist('Meta_Data')
else
[Partfile,Partpath] = uigetfile('*.mat','load the MetaData files');
                Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Channel_name);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Plot=0; % individual plotting
V1V=0; % individual tirals for Cross corr comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% User input

god=input('Which Site do you want to take a look at?     ');


Amp_interest = input('What amplitude?     ');


%% Choose the channesl that you want to look into
% Assuming Meta_Data is already loaded in your workspace
% Meta_Data is a cell array where Meta_Data{1, n}.Muscle_name contains the muscle names

% Get the total number of channels
num_channels = size(Meta_Data, 2);

% Initialize arrays to hold channel numbers and corresponding muscle names
channel_numbers = 1:num_channels;
muscle_names = cell(1, num_channels);

% Populate the muscle names array
for n = 1:num_channels
    muscle_names{n} = Meta_Data{1, n}.Muscle_name;
end

% Display channel numbers and their corresponding muscle names
disp('Channel Numbers and Corresponding Muscle Names:');
for n = 1:num_channels
    fprintf('Channel %d: %s\n', n, muscle_names{n});
end

% Prompt user to input selected channel numbers
prompt = 'Enter the channel numbers you want to select, separated by spaces: ';
selected_channels_input = input(prompt, 's');

% Convert the input string to an array of numbers
Selected_channels = str2num(selected_channels_input); %#ok<ST2NM>

% Validate the selected channels
if any(~ismember(Selected_channels, 1:num_channels))
    error('Selected channels contain invalid numbers. Please select numbers between 1 and %d.', num_channels);
end

% Display the selected channels
disp('Selected Channels:');
disp(Selected_channels);

% Display the corresponding muscle names for the selected channels
disp('Muscle Names for Selected Channels:');
for i = 1:length(Selected_channels)
    ch = Selected_channels(i);
    fprintf('Channel %d: %s\n', ch, muscle_names{ch});
end
% 
 


for k =1:1:length(Selected_channels)
   
    Channels{k} = Meta_Data{1,Selected_channels(k)};
end



%% Filter parameters
Wn = 45 / (Channels{1}.Snips_fs / 2);
[b, a] = butter(4, Wn); % Butterworth filter
numChannels=size(Channels,2);
% Initialize arrays to hold train data and info
Trains = cell(1, numChannels);
TrainsFiltered_Mean = cell(1, numChannels);
Train_info = cell(1, numChannels);

for k = 1:numChannels
    % Extract data for each channel
    [row, ~, ~] = find(Channels{k}.Guide_Matrix(:, 6) == Amp_interest);
    Sub_matrix_info = Channels{k}.Guide_Matrix(row, :);
    Sub_matrix_Trains = Channels{k}.Guide_trains(row, :);
    
    % Finding the correct sites
    [ro, ~, ~] = find(Sub_matrix_info(:, 1) == god);
    Train_info{k} = Sub_matrix_info(ro, :);
    Trains{k} = Sub_matrix_Trains(ro, :);
    
   
end

time = -500:1000/Channels{1}.Snips_fs:1000;

%% Plotting
figure;
for k = 1:numChannels
    subplot(2, 1, 1);
    plot(time, mean(Trains{k}), 'DisplayName', sprintf(' %s', Channels{k}.Muscle_name));
    hold on;
    subplot(2, 1, 2);
    plot(time, filtfilt(b,a,mean(Trains{k},1)), 'DisplayName', sprintf('Filt- %s', Channels{k}.Muscle_name));
     hold on;

   title(sprintf('site%d amp%g-Raw%g', god, Amp_interest,Unblank_raw));
    
end

       for hh=1:1:2
        subplot(2, 1, hh);
         xline(0, "LineStyle", "--");
         xline(500, "LineStyle", "--");
          grid on;
         legend;
        title(sprintf('site%d amp%g-Raw%g', god, Amp_interest,Unblank_raw));
        end



if Plot==1


% Individual trial plots
for i = 1:size(Train_info{1}, 1)
   J=figure;
     
        for k = 1:numChannels
            subplot(2, 1, 1);
            plot(time, Trains{k}(i, :), 'DisplayName', sprintf('%s', Channels{k}.Muscle_name));
            hold on;
            
            subplot(2, 1, 2);       
            plot(time,  filtfilt(b,a,Trains{k}(i, :)), 'DisplayName', sprintf('%s', Channels{k}.Muscle_name));
            hold on          
        end
        for hh=1:1:2
        subplot(2, 1, hh);
         xline(0, "LineStyle", "--", 'DisplayName','Onset');
         xline(500, "LineStyle", "--", 'DisplayName','OFFset');
          grid on;
         legend;
          title(sprintf('Trial %d-Raw%g', i,Unblank_raw));
        end

        %uiwait(J)
end

end
% 


%% compare any two signals
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I would choose Median over mean all the time for unfiltered signal
% for filtered signal, I would choose mean instead

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% with filter:
plotAlignedSignals(Channels, Selected_channels,Trains,time,'mean','filter',[])
plotAlignedSignals(Channels, Selected_channels,Trains,time,'median','filter',[])


%no filter:
plotAlignedSignals(Channels, Selected_channels,Trains,time,'mean','!',[])
plotAlignedSignals(Channels, Selected_channels,Trains,time,'median','!',[])

%% Trial by trial correlation

if V1V==1
num_channels2 = size(Meta_Data, 2);

% Extract muscle names and channel numbers
muscle_names2 = arrayfun(@(n) Meta_Data{1, n}.Muscle_name, 1:num_channels2, 'UniformOutput', false);

% Display channel numbers and their corresponding muscle names
fprintf('Channel Numbers and Corresponding Muscle Names:\n');

for n = 1:num_channels
    fprintf('Channel %d: %s\n', n, muscle_names2{n});
end

% Prompt user to input selected channel numbers
selected_channels_input2 = input('Enter exactly two channel numbers, separated by spaces: ', 's');
Selected_channels2 = str2num(selected_channels_input2);
% Check if exactly two channels are selected
if length(Selected_channels2) ~= 2 || any(~ismember(Selected_channels2, 1:num_channels))
    error('Please select exactly two valid channels between 1 and %d.', num_channels);
end

% Display selected channels and their muscle names
disp('Selected Channels:');
disp(Selected_channels2);
disp('Muscle Names for Selected Channels:');
for ch = Selected_channels2
    fprintf('Channel %d: %s\n', ch, muscle_names{ch});
end

% Extract data for the selected channels
Channels2 = Meta_Data(1, Selected_channels2);

Trains2 = cell(1, size(Channels2,2));

for k = 1:size(Channels2,2)
    % Extract data for each channel
    [row, ~, ~] = find(Channels2{k}.Guide_Matrix(:, 6) == Amp_interest);
    Sub_matrix_info2 = Channels2{k}.Guide_Matrix(row, :);
    Sub_matrix_Trains2 = Channels2{k}.Guide_trains(row, :);
    
    % Finding the correct sites
    [ro, ~, ~] = find(Sub_matrix_info2(:, 1) == god);
    Train_info{k} = Sub_matrix_info2(ro, :);
    Trains2{k} = Sub_matrix_Trains2(ro, :);
    
end



%%%
for i=1:size(Trains2{1,1},1)
plotAlignedSignals(Channels2, Selected_channels2,Trains2,time,'indiv','filter',i)

end

end


function plotAlignedSignals(Channels, Selected_channels,Trains,time,type,filter,Train_count)
    % Ensure selected_channels is a column vector
    num_selected = length(Selected_channels);
    Wn = 20 / (Channels{1}.Snips_fs / 2);
    [b, a] = butter(4, Wn); % Butterworth filter
    % Initialize figure for the matrix of plots
    Fu=figure;
    
    % Loop over each pair of selected channels
    for i = 1:num_selected
        for j = 1:num_selected
            % Select the current pair of channels
            ch1 = i;
            ch2 = j;
            
            % Extract signals'
            if strcmp(type, 'mean')
            signal1 =mean(Trains{i});
            signal2 =mean(Trains{j});
            elseif strcmp(type, 'median')
            signal1 =median(Trains{i});
            signal2 =median(Trains{j});
            elseif strcmp(type, 'indiv')
                 signal1 =Trains{i}(Train_count,:);
                 signal2 =Trains{j}(Train_count,:);
            end

            if strcmp(filter, 'filter')
               signal1= filtfilt(b,a,signal1);
               signal2= filtfilt(b,a,signal2);
            else
            end
%            
            % Calculate cross-correlation
            [xcorr_vals, lags] = xcorr(signal1, signal2, 'coeff');
            
            % Find the lag that maximizes the cross-correlation
            [~, idx] = max(xcorr_vals);
            lag = lags(idx);
            
            % Align signals based on the lag
%             if lag > 0
%                 aligned_signal2 = [nan(abs(lag),1); signal2(1:end-lag)];
%                 aligned_signal1 = signal1;
%             else
%                 aligned_signal1 = [nan(abs(lag),1); signal1(1:end-lag)];
%                 aligned_signal2 = signal2;
%             end
            
            % Create subplot for the current pair of channels
            subplot(num_selected, num_selected, (i-1)*num_selected + j);
            plot(time, circshift(signal2, round(lag)),'r', 'DisplayName', sprintf('Shifted %s', Channels{1,j}.Muscle_name));
            hold on;
            plot(time, signal1,'k', 'DisplayName', Channels{1,i}.Muscle_name);

           title(sprintf('Lag:%g %s vs %s %s Trial%g-Raw%g', 1000*(lag/Channels{1}.Snips_fs), Channels{1,i}.Muscle_name, Channels{1,j}.Muscle_name, type,Train_count),Unblank_raw);
%             xlabel('Sample Index');
%             ylabel('Signal Value');
            legend;
            grid on;

            
        end
        
    end
    if strcmp(type, 'indiv')
                uiwait(Fu)
    end
end



