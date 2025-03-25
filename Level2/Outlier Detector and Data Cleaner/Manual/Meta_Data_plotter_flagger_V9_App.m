%% load the proper channel
clear all

global  Unblank_raw Selected_Channel Snips_fs Guide_Matrix Threh_Dropdown_Value Arming Amp_interest  god Edit Mean_Cell     Median_Cell  Trains Train_info  Selected_Train_number_cell  Selected_Train_Matrix 


[Partfile,Partpath] = uigetfile('*.mat','load the MetaData files');
                Load_Channel_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Channel_name);




%% load and ask if you are editing or fresh start
 Edit=input(' Start fresh (Enter 0) OR  Edit an existing file (enter any key)    ');
if Edit==0
    Status='First_run';
    mkdir CleanedData
   Selected_Train_Matrix=nan(size(Meta_Data{1, 1}.Guide_Matrix,1),16);
   
else
    Status='Edit';
[Partfile,Partpath] = uigetfile('*.mat','Load the previous cleaned up data and continue where you left off');
                Load_Train_name=sprintf('%s\%s',Partpath,Partfile);
                load(Load_Train_name);
end



%% Preset Data

god=1;
% unique(Guide_Matrix(:,6))
Amp_interest=1;
Threh_Dropdown_Value=1;
Selected_Channel=1;
Guide_Matrix  = Meta_Data{1, Selected_Channel}.Guide_Matrix;
   


 createPopupDialog()
 



%% Staring a site Loop
while 1>0  % this is an infinite loop
% Finding the correct amplitudes
      Guide_trains  = Meta_Data{1, Selected_Channel}.Guide_trains ; 
      Guide_Matrix  = Meta_Data{1, Selected_Channel}.Guide_Matrix;
      Muscle_name   = Meta_Data{1, Selected_Channel}.Muscle_name;
              Case  = Meta_Data{1, Selected_Channel}.Case;
    Channel_Number  = Meta_Data{1, Selected_Channel}.Channel_Number;
          Snips_fs  = Meta_Data{1, Selected_Channel}.Snips_fs;

     
     [row,col, val] =find (Guide_Matrix(:,6)==Amp_interest);

     clear Sub_matrix_Trains Sub_matrix_info
     Sub_matrix_info= Guide_Matrix(row,:);
     Sub_matrix_Trains=Guide_trains(row,:);


     


    
    clear Trains   Train_info
    % finding the correct sites


    [ro,co,v]=find (Sub_matrix_info(:,1)==god);


if isempty(ro)
    

 createPopupDialog()
 
else



    Train_info=Sub_matrix_info(ro,:);
    Trains=Sub_matrix_Trains(ro,:);
     
    
    
    % normalize and plot in stacks
             


                 
    F= figure('units','pixels','outerposition', [0 0 816 960],'InnerPosition',[0 0 816 960])

            posFig = getpixelposition(F);
             hFig   = posFig(4);   % height of the figure
             Working_Step =0.85*hFig/size(Train_info,1);


             if Working_Step < 20
                 Steps=22;
                 Multirow= floor((22*size(Train_info,1))/(0.8*hFig));
             elseif Working_Step>=20
                 Steps=floor(0.85*Working_Step);
                Multirow=1;
             end
          

                clear stacked chkbox
                for i=1:1:size(Trains,1)
        
                        signal = Trains(i, :);

                        stacked=normalize(signal,"range")+ i;

                        pl= plot(stacked, 'DisplayName', sprintf('Tot%g-#%d',Train_info(i,3), Train_info(i,8)));
                        hold on
                        text(length(stacked)+2,i,sprintf('#%d-T%g',Train_info(i,8),Train_info(i,3)),"FontSize",8)
                        hold on   

                        
                           
                            if (0.1*hFig+i*Steps) >= (Working_Step*size(Trains,1))  &&  (0.1*hFig+i*Steps) < 2* (Working_Step*size(Trains,1))
                                 X_shift=70;
                                 Colomn_shift=1;       
                                 Y_loc= 0.1*hFig+(i-Colomn_shift*Current_i)*Steps;
                            elseif (0.1*hFig+i*Steps) >= 2* (Working_Step*size(Trains,1))
                                X_shift=2*70;
                                 Colomn_shift=2;  
                                 Y_loc= 0.1*hFig+(i-Colomn_shift*Current_i)*Steps;
                             elseif (0.1*hFig+i*Steps) >= 3* (Working_Step*size(Trains,1))
                                 X_shift=3*70;
                                 Colomn_shift=3; 
                                 Y_loc= 0.1*hFig+(i-Colomn_shift*Current_i)*Steps;
                            else 
                                X_shift=0;
                                Current_i=i;
                                Y_loc= 0.1*hFig+(i)*Steps;
                            end
                                
                            
                               
                       
                        chkbox(i) = uicontrol('Style','checkbox','String',sprintf('%g',i),            ...
                          'Value',0,'Position',[20+X_shift, Y_loc, 40, 20], ...
                          'Callback',{@checkBoxCallback,i});
                end
               
                  % Add a Saved indicator
          
                    Indicator=uicontrol('Style', 'text',  ...
                        'Position', [500, 0.05*hFig, 200, 30],'FontSize', 10, 'ForegroundColor','r');
                      Indicator.String=  {'Not Saved!'};


                % Here, I read the data from the files and if there is
                % anything, I want the checkboxes to be ticked. I do it
                % before select all becasue I do not want to deal with it.
                if exist('Selected_Train_Matrix')

                Preactivated_boxes=  ~isnan(Selected_Train_Matrix(Train_info(:,3),Channel_Number));
                [roo,Coll,Vals]= find(Preactivated_boxes==1);

                              
                    
                    % Set the value of each checkbox based on the "Select All" checkbox
                    for H = roo
                        if H <= length(chkbox)
                            set(chkbox(H), 'Value', 1);
                        end
                    end
                
                     Condition=sprintf('Modified_site: Included: %g out of %g',size(roo,1),size(Preactivated_boxes,1))
                else

                    Condition='Not_Modified'
                end




                  % Add "Select All" checkbox at the top
                    selectAllCheckbox = uicontrol('Style', 'checkbox', 'String', 'Select All', ...
                        'Position', [20, 0.05*hFig, 100, 20], ...
                        'Callback', {@selectAllCallback, size(Trains, 1)});
                  
                   % Add a button to compute and plot the mean of selected rows
                    uicontrol('Style', 'pushbutton', 'String', 'Plot Mean', ...
                        'Position', [200, 0.05*hFig, 100, 30], ...
                        'Callback', @(src, event) plotMeanCallback(Trains,Train_info,Case,Muscle_name));
                   
                    % Add a button to save the changes
                    uicontrol('Style', 'pushbutton', 'String', 'Save changes', ...
                        'Position', [400, 0.05*hFig, 100, 30], ...
                        'Callback', @(src, event) SaveChanges(Amp_interest,Channel_Number,god,Trains,Train_info,Case,Indicator));
                    
                        
                       % Add a Drop-down menu for the site number
                        Site_set=num2str(unique(Guide_Matrix(:,1)));
                        Selected=uicontrol('Style', 'popupmenu', 'Tag','Sites','String', Site_set, ...
                        'Position', [320, 0.05*hFig+5, 40, 30],'Callback',@(src, event) selection());
                          set(Selected,'Value',god)
                        % Add a Drop-down menu for the Threshold
                        Thresh_set=num2str(unique(Guide_Matrix(:,6)));
                        Threshsel=uicontrol('Style', 'popupmenu','Tag','Ampil', 'String', Thresh_set, ...
                        'Position', [110, 0.05*hFig+15, 40, 30],'Callback',@(src, event) selection_thresh());
                        set(Threshsel,'Value',find(str2num(Thresh_set)==Amp_interest))
                        % Add a Arm/Disarm switch                  
                        uicontrol('Style', 'togglebutton', 'String' , 'DisArmed', ...
                        'Position', [800, 0.05*hFig, 50,20 ],'BackgroundColor','w','Callback',@(src, event) Armer());
                       
                         % Add a Drop-down menu for the Channels
                        Channels=num2cell(1:1:size(Meta_Data,2));
                        Channels_sel=uicontrol('Style', 'popupmenu','Tag','chan', 'String', Channels, ...
                        'Position', [110, 0.05*hFig-15, 40, 30],'Callback',@(src, event) selection_Channel());
                        set(Channels_sel,'Value',Selected_Channel)




                 xline(Train_info(1,9),"LineStyle","--")
                 xline(Train_info(1,10),"LineStyle","--")

                 title(sprintf('(%s) STACKED Ch%g %s Real-Site%g Amp%g \n (Select non-outlier trials to keep) \n  Verdict: %g',Case,Channel_Number, Muscle_name,god,Train_info(1,6),Train_info(1,12)*Train_info(1,13)));
                 axis off
                 xlabel('Index');
                 ylabel('Signal Amplitude');

                 set(gcf, 'Position', get(0, 'Screensize'))% make it full screen

                uiwait(F)
                 
%                  input('must pause here')

end

end % end of god (for site numbers)


function Armer()
global Arming
Toggle = findobj(gcf, 'Style','togglebutton' );
         if Toggle.Value == 0
             Toggle.BackgroundColor= [1,1,1];
             Toggle.String='DisArmed';
              Arming=0;  
         elseif Toggle.Value == 1
              Toggle.BackgroundColor= [1,0,0];
              Toggle.String='Armed!';
              Arming=1;  
         end
             
end



function selection_thresh()
global Amp_interest Threh_Dropdown_Value
    Threshsel=findobj(gcf, 'Style','popupmenu','Tag','Ampil' );
    options=str2num(Threshsel.String);
    Amp_interest =options(Threshsel.Value);
    Threh_Dropdown_Value=Threshsel.Value;
    disp(Amp_interest)
    close all

end





function selection_Channel()
global   Selected_Channel
    Channel_sel=findobj(gcf, 'Style','popupmenu','Tag','chan' );
    %options=str2num(Channel_sel.String);
    Selected_Channel =Channel_sel.Value;
%     Chan_Dropdown_Value=Channel_sel.Value;
    disp(Selected_Channel)
    close all

end





function selection()
global god 
    DropDown = findobj(gcf, 'Style','popupmenu','Tag','Sites' );
    val = DropDown.Value;
    
    close all
    god=val;
    
    disp(val);
end



function checkBoxCallback(source, data, index)
     fprintf('Check box %s is %d\n',get(source,'String'),get(source,'Value'));
end


% Callback function for "Select All" checkbox
function selectAllCallback(source, ~, numCheckboxes)
    % Get the value of the "Select All" checkbox
    selectAllValue = get(source, 'Value');
    disp(selectAllValue)
    % Find all checkboxes in the current figure
    chkboxes = findobj(gcf, 'Style', 'checkbox');
    
    % Set the value of each checkbox based on the "Select All" checkbox
    for i = 1:numCheckboxes+1
        if i <= length(chkboxes)+1
            set(chkboxes(i), 'Value', selectAllValue);
        end
    end
end


function plotMeanCallback(Train_set,Traininfo,Case,Muscle)
global Snips_fs
    chkboxes = findobj(gcf, 'Style', 'checkbox');
    clear Active_Numbers
    for i=1:1:size(chkboxes,1)
        if isempty(str2num(get(chkboxes(i),'String'))) || get(chkboxes(i),'Value')==0
        else
            Active_Numbers(i)=str2num(get(chkboxes(i),'String'));
        end
    end
    nonZeroValues = sort(nonzeros(Active_Numbers));

    if isempty(nonZeroValues)
        disp('error! You must select at least one box')
    else
        
        MEAN=mean(Train_set(nonZeroValues,:),1);
        Median=median(Train_set(nonZeroValues,:),1);
        Total_Mean=mean(Train_set,1);
        Total_Median=median(Train_set,1);
        
        figure
                for j=1:1:2
                    if j==1
                    signals=MEAN;
                       titles='mean';
                    old_signal=  Total_Mean; 
                    elseif j==2

                     signals= Median ;
                     titles='median';
                     old_signal=  Total_Median; 
                    end
                %filtration
%                                
                                Wn=45/(Snips_fs/2); %Wn=Butter_Cut_Off/(Snips_fs/2); % Filter parameters
                                [b,a] = butter(4,Wn); % Set as butterworth filter
                                Butter_signal= filtfilt(b,a,signals);

             
                        subplot(3,1,j)
                        
                        rectangle('Position',[Traininfo(1,9),-10e9,(Traininfo(1,10)-Traininfo(1,9)),20e9],'FaceColor', ([227 229 229]/256),'EdgeColor','none')
                        
                        hold on

                        plot(normalize(Butter_signal,"range"),'k','LineWidth',0.5,'DisplayName',titles)
                        hold on
                        plot(normalize(filtfilt(b,a,old_signal),"range"),'r','LineWidth',0.5,'DisplayName','Selected All')
                        hold on
%                         text(0,0,sprintf('(%s) Mode:%s, Site:%g; Ch%g uA:%g,butter=%g',Case,Mode,i,Channel_of_intrest,Amplitude_of_interest,Butter_Cut_Off), "color",'r' )
                       %%%%%%%
                       ylim([0 1.2])          % [mean(signals)-5*(std(signals)) 1.05*max(signals)]) %ylim([mean(MEAN)-5*(std(MEAN)) mean(MEAN)+5*(std(MEAN))])
                       %%%%%%%
                       axis off
                        set(gcf, 'color', 'w'); 
                        legend
                        
                        title(sprintf('(%s) Site:%g uA:%g Verdict:%guA Ch:%s %s \n Included: %g  Excluded: %g ',Case,Traininfo(1,1),Traininfo(1,6),Traininfo(1,12)*Traininfo(1,13),Muscle,titles,size(nonZeroValues,1),size(Traininfo,1)- size(nonZeroValues,1)));   %%%%inc:%g Rej:%g
                        
                        
                end
                       subplot(3,1,3)
                        rectangle('Position',[Traininfo(1,9),-10e9,(Traininfo(1,10)-Traininfo(1,9)),20e9],'FaceColor', ([227 229 229]/256),'EdgeColor','none')
                        
                        hold on

                       plot(normalize(filtfilt(b,a,MEAN),"range"),'DisplayName','Mean')
                       hold on 
                       plot(normalize(filtfilt(b,a,Median),"range"),'DisplayName','Median')
                       hold on
                        ylim([0 1.2])          % [mean(signals)-5*(std(signals)) 1.05*max(signals)]) %ylim([mean(MEAN)-5*(std(MEAN)) mean(MEAN)+5*(std(MEAN))])
                       %%%%%%%
                       axis off
                        set(gcf, 'color', 'w'); 
                       legend
                         title(sprintf('(%s) Site:%g uA:%g Verdict:%guA Ch:%s Mean vs Median \n Included: %g  Excluded: %g ',Case,Traininfo(1,1),Traininfo(1,6),Traininfo(1,12)*Traininfo(1,13),Muscle,size(nonZeroValues,1),size(Traininfo,1)- size(nonZeroValues,1)));   %%%%inc:%g Rej:%g
                       
    end
    
end



                

function SaveChanges(Amp,Chan,sitenum,Trainset,TrainInfo,CASE,Indicator)
    global Arming Edit Mean_Cell     Median_Cell    Selected_Train_number_cell   Selected_Train_Matrix Unblank_raw

    chkboxes = findobj(gcf, 'Style', 'checkbox');
    clear Active_Numbers
    for i=1:1:size(chkboxes,1)
        if isempty(str2num(get(chkboxes(i),'String'))) || get(chkboxes(i),'Value')==0
        else
            Active_Numbers(i)=str2num(get(chkboxes(i),'String'));
        end
    end
    nonZeroValues = sort(nonzeros(Active_Numbers));

    if isempty(nonZeroValues)
        disp('error! You must select at least one box')
    else
        

        Mean_Cell{Amp,Chan}(sitenum,:)= mean(Trainset(nonZeroValues,:),1);
        Median_Cell{Amp,Chan}(sitenum,:)= median(Trainset(nonZeroValues,:),1);
        if Edit==0
        else
            Selected_Train_number_cell{Amp,Chan}{sitenum,1}=[];
        end
            Selected_Train_number_cell{Amp,Chan}{sitenum,1}=TrainInfo(nonZeroValues,3); % this is absolute number

        
        Selected_Train_Matrix(TrainInfo(1,3):TrainInfo(end,3),Chan) =nan;
          
        % active number to real number conversion
        

        for j=nonZeroValues
            
        Selected_Train_Matrix(TrainInfo(j,3),Chan)= 1;
        end

        % MAKE SURE TO SAVE THE SYSTEM with version control
        if Arming==1
        Type='Manual_Clean';
        save(sprintf('CleanedData/%s_CleanedData_(%s).mat',CASE,datestr(now, 'dd-mmm-yyyy')),'Mean_Cell','Median_Cell','Selected_Train_number_cell','Selected_Train_Matrix','Type','Unblank_raw')
        
        Indicator.String=  {'Done! and Saved!'};
        Indicator.ForegroundColor=[0 1 0];
        end
       
    end


end


function createPopupDialog()
global god Amp_interest Guide_Matrix

    % Create a figure for the dialog
    fig = uifigure('Name', 'Popup Dialog Box', 'Position', [100, 100, 400, 300]);

    % Create a label for the text message
    lblMessage = uilabel(fig, 'Position', [50, 220, 300, 30], 'Text', sprintf('Site %g, Amp %g does not exist. \n Please select another threshold or site:',god,Amp_interest));
    
    % Create the drop-down menu for Site Number
    lblSite = uilabel(fig, 'Position', [50, 160, 100, 30], 'Text', 'Site Number:');
    ddSite = uidropdown(fig, 'Position', [160, 160, 150, 30], 'Items', cellstr(num2str(unique(Guide_Matrix(:,1))))  );
 
    % Create the drop-down menu for Amplitude
   
    lblAmplitude = uilabel(fig, 'Position', [50, 110, 100, 30], 'Text', 'Amplitude:');
    ddAmplitude = uidropdown(fig, 'Position', [160, 110, 150, 30], 'Items',cellstr(num2str(unique(Guide_Matrix(:,6)))) );
    
    
    % Create the Go button, centered
    btnGo = uibutton(fig, 'push', 'Text', 'Go', 'Position', [160, 40, 80, 30], 'ButtonPushedFcn', @(src, event) goButtonPressed());

    % Add callbacks for value changes
    uiwait(fig)
    function goButtonPressed()
        % Display selected values in the Command Window
           god=str2num(ddSite.Value);
           Amp_interest=str2num(ddAmplitude.Value);
        close(fig)
    end
end
