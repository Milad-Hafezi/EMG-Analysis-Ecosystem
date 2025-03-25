%% voronoi diagram generation

% SVG output is also added

% New colors added after undetrending:
% Color_map([5 10 20 40 80 160 320],:)=([13 19 61;38 34 98; 13 19 61;41 44 106; 16 115 181;54 186 220;163 222 249])/256;
% Changed colors for 20uA-320uA because undetrended verdict has no
% 5-10uA-positive sites. 5-10uA colors are unchanged. (JL edit)




% Ai to do:


            % * all the concaves and the curves must be cleaned up
            % *



% needs to be done:

            % (DONE)*axis off 
            % (DONE)*outline should be colored differnt than blue

            % (DONE) *any other cells that are not active to Grey      RGB:
            %  [0.909803921568627   0.905882352941176   0.898039215686275] 

      

            % (DONE)*darker shades of grey for not a fair test
            % (DONE)*put swtich for number labels and center dots


            % *Think about how boolan maps could be plotted at once so that
            %  it would be easier for the use to use them later on .

            % (DONE) Make a plot for how many muscles are activated per each site
            % with shades of green.


            % set 2 to dos:

            % - Shaded boolean based on Dylan's idea. make 7 shades or
            % colors and find the corresponding thresholds and build sites
            % that have lower than a limit threhsold



            %- (Done)Jaspreet's help for relative locations of the 




close all
% clear all
clearvars -except cells Extract Site_Row_Key Scale_info %[x1 y1; x2 y2; min(x1,x2)  pixelScaleDist+y1];

Case='ERROR'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Boolean_Logic=input('Do you want to run the booleans too?  0 for no and 1 for yes     ');  % 1 to generate the boolean plots 0 to skip them
Devel=0;
Number_logic=0;
Scale_lining=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Thresh_sequence  = [320 160 80 40 20 10 5];  % this is used as all the amps that will be plotted in booleans

%% Make Folders

mkdir Heatmaps_fig     
mkdir Heatmaps_eps 
mkdir Booleans
mkdir Comparison
mkdir Centroid_logs

    mkdir Heatmaps_SVG Booleans
    mkdir Heatmaps_SVG Heatmaps
    mkdir Heatmaps_SVG Greens
    mkdir Heatmaps_SVG Compare
    

    


%% Color pot

Grey_Back=[0.909803921568627   0.905882352941176   0.898039215686275];
Grey_lines=  [0.819607843137255   0.827450980392157   0.831372549019608];
Grey_Unfair= [0.498039215686275   0.498039215686275   0.498039215686275];

%% loading the Verdict

[file,path] = uigetfile('*.xlsx','Select the Verdict file for your case');
cd(path);
Verdict=importdata(file);

MUSCLES=Verdict.textdata.Sheet17(1:16,1)
Case=Verdict.textdata.Sheet17{1, 2} 
Color=Verdict.data.Sheet17(1:16, 2:4) 
Color=Color./256;

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
    Coactivation(:,Bin)=C{Bin}(:,3);
    CoThresh(:,Bin)=C{Bin}(:,4);
end

for Jin=1:1:size(Coactivation,1)
    NNZ(Jin)=nnz(Coactivation(Jin,:));
end
NNZ=NNZ';

%% Ploting the Vertices
     fig = figure;
     hold on;
     for id = 1:size(Extract,1)
        %plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
        p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
        p.EdgeColor=Grey_lines;
%        plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
       
%        text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
               
     end
     if Scale_lining==1
          line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
          line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
     end
    axis off
    axis equal
        
    xl = xlim; % getting the X and Y limits to put in for later plots
    yl= ylim;


%% Heat Maps
color1 = [0.5, 0.5, 1];      % Set colors according to your definition of
color2 = [0,   0,   0.5];    % "light blue" and "dark blue"

Color_map([5 10 20 40 80 160 320],:)=([13 19 61;38 34 98; 13 19 61;41 44 106; 16 115 181;54 186 220;163 222 249])/256;
% Changed colors for 20uA-320uA because undetrended verdict has no
% 5-10uA-positive sites. 5-10uA colors are unchanged. (JL edit)
for i=1:1:size(C,1)
  
    clear XX YY THR
    XX=[];
    PILE=[];

    Channel=C{i}(1,1);
    kiv=1;
    for k=1:1:size(C{i},1)
        clear Redo Redos
 
        Site=C{i}(k,2);
              
        %%%%%%%%%%%%%%%%%%%%%%% Detecting the correct threshold
        
        [col loc]=find (cell2mat(Extract(:, 7))==Site);
        if ~isempty(col)
            Row=col;
            Skip=0;
                      Redos= cell2mat(Extract(Row,9));
                      if ~Redos== 0 % if there are no redos (Redos=0)
                       PILE=[PILE,Redos];
                       clear temp_thresh
                                 for j=1:1:length(Redos)
                                           [co lo]=find(C{i}(:,2)==Redos(j));
                                           tmp_th=C{i}(co,4);
    
                                           if tmp_th > 0
                                            temp_thresh(j)=tmp_th;
                                           else
                                            temp_thresh(j)=nan;
                                             end
                                 end
                       Thresh = min(temp_thresh);
                      elseif Redos== 0
                      Thresh =  C{i}(k,4); 
                      end
                         
          
        elseif isempty(col)
                    
                      %  we want to see if this redo bach has been dealt with
                         % before
              if   ismember(Site,PILE)
                        Skip=1;  %Skip if it has been dealt with before
              elseif  ~ismember(Site,PILE)
                        Skip=0; % donot skip
                      for jj=1:1:length(Extract(:,9)) % find the row that this redo is located in
                              [coll locc]= find(cell2mat(Extract(jj,9))==Site);
                                    if length(coll) > 0 && length (coll) < 2
                                        Row=coll;
                                    end
                      end
%                    Threhsold finding   
                   Redos=cell2mat(Extract(Row,9));                   
                   PILE=[PILE,Redos];
                   clear temp_thresh
                             for j=1:1:length(Redos)
                                       [co lo]=find(C{i}(:,2)==Redos(j));
                                       tmp_th=C{i}(co,4);

                                   if tmp_th > 0
                                    temp_thresh(j)=tmp_th;
                                   else
                                    temp_thresh(j)=nan;
                                     end
                             end
                   Thresh = min(temp_thresh);


              end
        end
                        

            
        if Skip == 0

                YESNO=C{i}(k,3);
                            if ~(Redos)==0 % if there are no redos (Redos=0)
            
                            YESNO=sum(C{i}(Redos,3));
                            end

                if YESNO>=1
        %         color = color1 * (Thresh/320) + color2 * (1 - Thresh/320);
        %         delta(k,:)=color
                    if Thresh ==0
                        Thresh=5;
                    end

                    if Thresh==888
                        color=Grey_Unfair;
                    else
                    color=Color_map(Thresh,:);
                    end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
             %% filling the votecies
               
                taco=Extract{Row, 6}(:,1);
                baco=Extract{Row, 6}(:,2);

                p=fill(taco ,-baco,color);
                p.EdgeColor=Grey_lines;
                % loging the active sites and their threshold with they location
                  [r c]=find(cell2mat(Extract(:,7))==Site);  % finding center coordinates
                if Thresh==88
                else
                XX(kiv)=cell2mat(Extract(r,2));
                YY(kiv)=-cell2mat(Extract(r,3));
                THR(kiv)=1/Thresh; % this is the excitability
                kiv=kiv+1;
                end
                end
        end
    end

    
        if size(XX,1) ==0 % if the whole channel is 0  we do this to avoid error
            XX=0;
            YY=0;
            THR=1;
        end



        % calculating the weighted agerage
        Mean_X=mean(XX);
        Mean_Y=mean(YY);
        %Weighted mean
        WMean_X= sum(XX.*THR)/sum(THR);
        WMean_Y= sum(YY.*THR)/sum(THR);

        WDev_X=sqrt(var(XX,THR));
        WDev_Y=sqrt(var(YY,THR));

        % Ellipse equation and determination.
           T=atan(Mean_Y-WMean_Y)/(Mean_X-WMean_X);
            clear x1 x2 y
           b=std(YY);
           a= std(XX);
             scatter(WMean_X,WMean_Y,'filled','MarkerFaceColor',Color(i,:)) 
             scatter(Mean_X,Mean_Y,'d','MarkerFaceColor',Color(i,:))
           
            x1=Mean_X-std(XX):1:Mean_X;
            x2=Mean_X:1:Mean_X+std(XX)+1;
            x=[x1 x2];
            y=((b/a)*sqrt((a^2)-(x-Mean_X).^2));
            
            plot(x,y-abs(Mean_Y),'r')
            plot(x,-y+Mean_Y,'r')
            hold on

            MishMoosh(Channel,1)=a;
            MishMoosh(Channel,2)=b;
            MishMoosh(Channel,3)=Mean_X;
            MishMoosh(Channel,4)=Mean_Y; 
            MishMoosh(Channel,5)=std(XX);
            MishMoosh(Channel,6)=std(YY);

            b=WDev_Y;
            a= WDev_X;
            clear x1 x2 y

            x1=WMean_X-WDev_X:1:WMean_X;
            x2=WMean_X:1:WMean_X+WDev_X+1;
            x=[x1 x2];
            y=((b/a)*sqrt((a^2)-(x-WMean_X).^2));
            
            plot(x,y-abs(WMean_Y),'--g','LineWidth',2)
            plot(x,-y+WMean_Y,'--g','LineWidth',2)

            WMishMoosh(Channel,1)=a;
            WMishMoosh(Channel,2)=b;
            WMishMoosh(Channel,3)=WMean_X;
            WMishMoosh(Channel,4)=WMean_Y;
            WMishMoosh(Channel,5)=WDev_X;
            WMishMoosh(Channel,6)=WDev_Y;
%%%%%%%%
    if Number_logic==1
        for id = 1:size(Extract,1)
         text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
        end
    end
        title(sprintf('%s Ch%g, %s', Case,Channel,MUSCLES{i}));
        savefig(sprintf('Heatmaps_fig/HeatMap_Channel%g_%s.fig',Channel,Case))
        exportgraphics(gca,sprintf('Heatmaps_eps/HeatMap_Channel%g_%s.eps',Channel,Case)) %exporting to eps
               saveas(gca,sprintf('Heatmaps_SVG/Heatmaps/HeatMap_Channel%g_%s.svg',Channel,Case))  % saveas SVG

        close all
%% plot again
          fig = figure;
                 hold on;
                 for id = 1:size(Extract,1)
%                        plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                         p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                         p.EdgeColor=Grey_lines;
            %        plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                   
            %        text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
                           
                 end

                 if Scale_lining==1
                    line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
                    line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
                 end


        hold on;
        axis off
        grid off
        axis equal
       xlim(xl)
       ylim(yl)
         
end



%% Boolean map 2nd order (and)
if Boolean_Logic ==1
for i=1:1:size(C,1)



    Channel=C{i}(1,1) 
    Yes_col=C{i}(:,3); %getting the channel of interest
    Thresh_Col=C{i}(:,4); %getting their thresholds

    clear Match
    for k=1:1:size(C,1)

         




       Yes_col_2nd=C{k}(:,3);
       Thresh_Col_2nd=C{k}(:,4);

       Pre_rubric=Yes_col.*Yes_col_2nd;
       rubric=  max(Thresh_Col,Thresh_Col_2nd).*Pre_rubric;  %sqrt(Thresh_Col.*Thresh_Col_2nd);

           for in=Thresh_sequence
                
        
               [r c]=find(rubric <= in & rubric > 0);
               Match=C{i}(r,2);
               %colora=  [1.000000000000000  0.666666666666667   0]* ((length(Thresh_sequence)+1- in)/length(Thresh_sequence)) ;
               colora= Color_map(in,:);
                for J=1:1:size(Match,1)
                    
                         Site=Match(J);
                        
                        [r ,c]= find(Site_Row_Key(:,1) == Site);
                         Row=Site_Row_Key(r,2);
                
                            taco=Extract{Row, 6}(:,1);
                            baco=Extract{Row, 6}(:,2);
                            p=fill(taco ,-baco,colora);
                            p.EdgeColor=Grey_lines;
                end
           end


            if Number_logic==1
             for id = 1:size(Extract,1)
             text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
             end
            end
         title(sprintf('Bool2 %s Ch%g VS Ch%g, %s vs %s', Case,Channel,k,MUSCLES{Channel},MUSCLES{k}));
%         savefig(sprintf('Booleans/%s Boolean Map Ch%g VS Ch%g.fig',Case,Channel,k))
        exportgraphics(gca,sprintf('Booleans/%s Boolean Map Ch%g VS Ch%g.eps',Case,Channel,k)) %exporting to eps
                saveas(gca,sprintf('Heatmaps_SVG/Booleans/%s Boolean Map Ch%g VS Ch%g.svg',Case,Channel,k)) % SVG
        close all
                %% plot again
                  fig = figure;
                         hold on;
                         for id = 1:size(Extract,1)
%                               plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                                p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                                 p.EdgeColor=Grey_lines;
                    %        plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                           
                    %        text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
                                   
                         end
                         if Scale_lining==1
                           line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
                           line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
                         end
        
                hold on;
                axis off
                grid off
                axis equal
                xlim(xl)
                ylim(yl)
    end
        
end

end

%% OR map  OR plot analysis
Smaple_Thresh= C{1}(:,4);
Valid_Threhsh=Smaple_Thresh < 888;
Sum = sum(Coactivation,2).*Valid_Threhsh;

for Chen=Thresh_sequence
        clear S
        for ink=1:1:size(CoThresh,1)
           
           [R, CC]= find(CoThresh(ink,:)<=Chen & CoThresh(ink,:) > 0);
           S(ink)= nnz(CC);
                  
        end
        S=S';


        for J=1:1:size(Extract,1)
               Site_of_in= Extract{J,7};
               Redocheck=Extract{J,9};
    
               if nnz(Redocheck) == 0
                   % no Redos
                   MuscleNum=S(Site_of_in);
                   Site_to_plot=Site_of_in;
                elseif nnz(Redocheck)>0
                   % we have Redos
                   [M,I]=max(S(Redocheck));
                   Site_to_plot=Redocheck(I);
                   MuscleNum=M;
               else
                   % Error
                   error('88888888888888888888888888')
               end

    
                taco=Extract{J, 6}(:,1);
                baco=Extract{J, 6}(:,2);

            if MuscleNum>0 
                colo= [0.501960784313725   1.00   0] *  ((16- MuscleNum)/16) ; % generating shades of green
                p=fill(taco ,-baco,colo);
                p.EdgeColor=Grey_lines;

                 text((Extract{J,2}-1),-(Extract{J,3}-1),num2str(MuscleNum),"FontSize",9,'Color','r')
            end
        end



        
        if Number_logic==1
         for id = 1:size(Extract,1)
         text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
         end
        end

    % Finding the redo sites with most muscles active
            
        title(sprintf('Green Plot %s , Amplitude (%g) uA', Case,Chen));
        savefig(sprintf('Comparison/%s Green_plot Amp %g.fig',Case,Chen))
        exportgraphics(gca,sprintf('Comparison/%s Green_plot Amp %g.eps',Case,Chen)) %exporting to eps
                saveas(gca,sprintf('Heatmaps_SVG/Greens/%s Green_plot Amp %g.svg',Case,Chen)) %  SVG


close all
 %% plot again
                  fig = figure;
                         hold on;
                         for id = 1:size(Extract,1)
%                               plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                                p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                                 p.EdgeColor=Grey_lines;
                    %        plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                           
                    %        text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
                                   
                         end
                         if Scale_lining==1
                           line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
                           line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
                         end
        
                hold on;
                axis off
                grid off
                axis equal
                xlim(xl)
                ylim(yl)




end
        


  %% Comparison Plot for mean
        
        %%          plot again
                  fig = figure;
                         hold on;
                         for id = 1:size(Extract,1)
                               plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                                 p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                                 p.EdgeColor=Grey_lines;
                    %        plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                           
                    %        text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
                                   
                         end

                         if Scale_lining==1
                            line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
                            line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
                         end
        
        
                    hold on;
                    axis off
                    grid off
                    axis equal
                    xlim(xl)
                    ylim(yl)
                    %%
            
        
        %          newcolors = [0.83 0.14 0.14; 1.00 0.54 0.00; 0.47 0.25 0.80; 0.25 0.80 0.54];
                  PP=[];
                  LG='';
                  bictor=1;
                for bic=1:1:size(MishMoosh,1)
                    clear x1 x2 x
                     a=MishMoosh(bic,1);
                     b= MishMoosh(bic,2);
                     Mean_X=MishMoosh(bic,3);
                     Mean_Y=MishMoosh(bic,4);
                     
        %              scatter(WMean_X,WMean_Y,'filled') 
                     
                    x1=(Mean_X-MishMoosh(bic,5)):1:Mean_X;
                    x2=Mean_X:1:(Mean_X+MishMoosh(bic,5)+1);
                    x=[x1 x2];
                    y=((b/a)*sqrt((a^2)-(x-Mean_X).^2));
                             
        %                     colororder(newcolors)
                    if Mean_X ==0 &&  Mean_Y==0
                        %do nothing
                    else
                        scatter(Mean_X,Mean_Y,'d','MarkerFaceColor',Color(bic,:))
                        p= plot(x,y-abs(Mean_Y),x,(-y+Mean_Y),'LineWidth',2,'Color',Color(bic,:));
                        text(Mean_X,Mean_Y, MUSCLES{bic}, 'FontSize', 8);
                        
                        PP=[PP,p];
                        LG{bictor}=sprintf('Ch%g , %s',bic,MUSCLES{bic})  %sprintf('Ch%g',bic);
                        bictor=bictor+1;
                        
                    end
                end
                   legend(PP(1,:),LG, 'Location','southoutside','Orientation','horizontal','NumColumns',4) 
                   title(sprintf('MishMoosh %s Activity Zones', Case));
                   savefig(sprintf('Comparison/%s Mishmoosh plot.fig',Case))
                   exportgraphics(gca,sprintf('Comparison/%s Mishmoosh plot.eps',Case)) %exporting to eps
                           saveas(gca,sprintf('Heatmaps_SVG/Compare/%s Mishmoosh plot.svg',Case)) %exporting to eps
                close all


    %% Comparison Plot for Weighted mean
   %%          plot again
                  fig = figure;
                         hold on;
                         for id = 1:size(Extract,1)
%                                plot([Extract{id,6}(:,1);Extract{id,6}(1,1)],-[Extract{id,6}(:,2);Extract{id,6}(1,2)],'color', Grey_lines, 'LineWidth', 1)
                                 p = fill(Extract{id, 6}(:,1) ,-Extract{id, 6}(:,2),Grey_Back);
                                 p.EdgeColor=Grey_lines;
                    %        plot(Extract{id,2},-Extract{id,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
                           
                    %        text((Extract{id,2}+1),-(Extract{id,3}+1),Extract{id,1},"FontSize",6) 
                                   
                         end
                         if Scale_lining==1
                            line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
                            line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
                         end
        
        
                    hold on;
                    axis off
                    grid off
                    axis equal
                    xlim(xl)
                    ylim(yl)
                    %%
%          newcolors = [0.83 0.14 0.14;0.83 0.14 0.14; 1.00 0.54 0.00; 1.00 0.54 0.00; 0.47 0.25 0.80; 0.47 0.25 0.80; 0.25 0.80 0.54;0.25 0.80 0.54];
         PP=[];
         LG='';
         bictor=1;
        for bic=1:1:size(WMishMoosh,1)
            clear x1 x2 x
             a=WMishMoosh(bic,1);
             b= WMishMoosh(bic,2);
             WMean_X=WMishMoosh(bic,3);
             WMean_Y=WMishMoosh(bic,4);
             
                   
                    x1=WMean_X-WMishMoosh(bic,5):1:WMean_X;
                    x2=WMean_X:1:WMean_X+WMishMoosh(bic,5)+1;
                    x=[x1 x2];
                    y=((b/a)*sqrt((a^2)-(x-WMean_X).^2));
                            
        %                     colororder(newcolors)
                    
        %             legend('Orientation','horizontal','Location','southoutside','NumColumns',16)
                    if WMean_X ==0 &&  WMean_Y==0
%                          %do nothing
                     else
                           scatter(WMean_X,WMean_Y,'filled','MarkerFaceColor',Color(bic,:)) 
                     
                         p=plot(x,y-abs(WMean_Y), x,(-y+WMean_Y),'LineWidth',2,'Color',Color(bic,:))
                            
                          text(WMean_X,WMean_Y, MUSCLES{bic}, 'FontSize', 8);
                    
                    PP=[PP,p];
                    LG{bictor}=sprintf('Ch%g , %s',bic,MUSCLES{bic})  %sprintf('Ch%g',bic);
                    bictor=bictor+1;
                    end
        %             plot(x,y-abs(WMean_Y),'LineWidth',2,'DisplayName',sprintf('Ch%g',bic))
        %             hold on
        %             plot(x,-y+WMean_Y)
        %             hold on
                   
                    
        end

           legend(PP(1,:),LG, 'Location','southoutside','Orientation','horizontal','NumColumns',4) % puuting on the legends
           title(sprintf('Weighted MishMoosh %s Activity Zones', Case));
           savefig(sprintf('Comparison/%s Weighted Mishmoosh plot.fig',Case))
           exportgraphics(gca,sprintf('Comparison/%s Weighted Mishmoosh plot.eps',Case)) %exporting to eps
                   saveas(gca,sprintf('Heatmaps_SVG/Compare/%s Weighted Mishmoosh plot.svg',Case)) %exporting to svg
        close all


        %% Save the location of elipses for weighted Centroids

        filename = sprintf('Centroid_logs/Weighted Centroid log Case %s.mat', Case);              
        Savable=WMishMoosh(:,1:4);  % cutting the first part of Wmishmoosh
        save(filename,"Savable")

        filename2 = sprintf('Centroid_logs/Weighted Centroid log Case %s.xlsx', Case);
        Savable_names= [num2cell(Savable) MUSCLES]
        sheet = 1;
        col_header={'WDev_X','WDev_Y','WMean_X','WMean_Y','Muscle name','Channel number'};
        Channel_List=num2cell([1:16]');
        output_matrix=[col_header; [Savable_names,Channel_List]]; 
        xlswrite(filename2,output_matrix,sheet)


        %% Save the location of elipses for Normal Centroids
  
       

        filename = sprintf('Centroid_logs/Normal Centroid log Case %s.mat', Case);              
        Savable=MishMoosh(:,1:4);  % cutting the first part of Wmishmoosh
        save(filename,"Savable")

        filename2 = sprintf('Centroid_logs/Normal Centroid log Case %s.xlsx', Case);
        Savable_names= [num2cell(Savable) MUSCLES]
        sheet = 1;
        col_header={'Dev_X','Dev_Y','Mean_X','Mean_Y','Muscle name','Channel number'};
        Channel_List=num2cell([1:16]');
        output_matrix=[col_header; [Savable_names,Channel_List]]; 
        xlswrite(filename2,output_matrix,sheet)

