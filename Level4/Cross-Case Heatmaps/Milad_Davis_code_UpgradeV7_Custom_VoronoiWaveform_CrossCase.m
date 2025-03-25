
%Stuff to do: 
            % 1-  (DONE) redo sites must be able to be seperated ( '55/70/90' must be considered as site 55 and site 70 and 90 as REDO) 
            % 2- (DONE)         We should be able to functionally plot all of the
                                % voronoies and color them as needed

            % 3- (DONE) Patching or filling the vortecies

            % 5- (DONE) We might be able to use the automatic center detector? In
            % that case, I have to double check the code

            %6- (Done) Figure out if you need to rescale





 [file,path] = uigetfile('*.svg','Select the SVG File');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Scaling_On_Off= 1;  % 0 for scaling "1 to 1" and 1 for autoscaling based on the scale line

scaleLength = 1000;% microns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(path); %file path is not updated
% fid = fopen('23-63 Voronoi_v3 (border cropped)MH Edit copy.svg')
fid = fopen(file);
textFile = textscan( fid, '%s', 'Delimiter', '\r' );

%% scaling
% get the scale bar coordinates
if Scaling_On_Off ==1  
        
        lineLines = find(contains(textFile{1,1},'<g id="ruler">'))+1; % finds the next line after <g id="ruler"> 
        % You should consistently use a the name "ruler" for you scale line title
        if size(lineLines)==1
            chunk = textFile{1,1}{lineLines,1};
            x1 = str2num(cell2mat(extractBetween(chunk,'x1="','" y1=')));
            y1 = str2num(cell2mat(extractBetween(chunk,'y1="','" x2=')));
            x2 = str2num(cell2mat(extractBetween(chunk,'x2="','" y2=')));
            y2 = str2num(cell2mat(extractBetween(chunk,'y2="','"/>')));
        elseif size(lineLines,1)>1
            try
                scaleLinesStart = find(contains(textFile{1,1},'<g id="scale'));
            catch
                scaleLinesStart = find(contains(textFile{1,1},'<g id="Scale'));
            end
            nextG = find(contains(textFile{1,1}(scaleLinesStart:end,:),'</g>'));
            scaleLinesEnd = scaleLinesStart+nextG(1)-1;
            chunk = strjoin(textFile{1,1}(scaleLinesStart:scaleLinesEnd,:));
            x1 = str2num(cell2mat(extractBetween(chunk,'x1="','" y1=')));
            y1 = str2num(cell2mat(extractBetween(chunk,'y1="','" x2=')));
            x2 = str2num(cell2mat(extractBetween(chunk,'x2="','" y2=')));
            y2 = str2num(cell2mat(extractBetween(chunk,'y2="','"/>')));
            
        elseif size(lineLines,1)<1
            error('No lines! Make sure your image contains a scale bar')
        end
        % get the pixels-microns conversion factor
        try
            pixelScaleDist = sqrt((x2-x1)^2+(y2-y1)^2); % length of scale bar in pixels
        catch
            error('Too many lines, and no labeled scale bar! Make sure either the only line in your image is a scale bar, or you label the layer with "scale"')
        end
        
        resolution = scaleLength/pixelScaleDist; % microns per unit in original data
        
        

        Scale_info=[x1 -y1; x2 -y2; min(x1,x2)  -y1+pixelScaleDist]; % saving the point of the scale line for furthur use NOTE: the third line is the Perpendicular line
        Scale_info= resolution* Scale_info;
        
elseif Scaling_On_Off==0

        resolution=1;
end
  

%% Labels or sites
% Prep: make a matrix of label x and y coordinates, and ID
% find all lines that start with "<text"
textLines = find(contains(textFile{1,1},'<text'));
circles = find(contains(textFile{1,1},'<circle'));

%look for circles for site markers first


if ~isempty(textLines)
    for i = 1:numel(textLines)
        chunk = textFile{1,1}{(i),1};
        if ~isempty(circles)
            for i = 1:numel(circles)
                chunk = textFile{1,1}{circles(i),1};
                xC = str2num(cell2mat(extractBetween(chunk,'cx="','" cy="')));
                yC = str2num(cell2mat(extractBetween(chunk,'" cy="','" r="')));
                xC = xC*resolution;
                yC = yC*resolution;
                possibleSites{i,1} = xC;
                possibleSites{i,2} = yC;
                possibleSites{i,3} = NaN;
                easyPlotSites(i,1) = xC;
                easyPlotSites(i,2) = yC;
            end
        end




        if ~isempty(textLines)
            % find the center of each polygon further down
            if size(textLines,1)>=numel(circles)
                possibleLabels = {NaN(numel(textLines),3)};
                for i = 1:numel(textLines)
                    chunk = textFile{1,1}{textLines(i),1};
                    
                    % check if the line ends with '>'
                    endchar = strfind(chunk,'>');
                    
                    % if not, then grab the next line down and concatenate
                    if isempty(endchar)
                        chunk2 = textFile{1,1}{textLines(i)+1,1};
                        chunk = [chunk chunk2];
                    end
                    
                    % x-coordinate is after "matrix(1 0 0 1" and before " "
                    if contains(chunk, 'matrix(1 0 0 1 ' ) == 1
                        x = cell2mat(extractBetween(chunk,'matrix(1 0 0 1 ',' '));
                        % y-coordinate is after x-coordinate, after " ", and before ")"
                        y = cell2mat(extractBetween(chunk,[x ' '],')'));
                    else
                        chunk1 = cell2mat(extractBetween(chunk,'matrix(',')'));
                        spaceLoc = isspace(chunk1);
                        x = find(spaceLoc == 1);
                        x = x(:, 4);
                        x = cell2mat(extractBetween(chunk1,x+1,' '));
                        
                        % y-coordinate is after x-coordinate, after " ", and before ")"
                        y = cell2mat(extractBetween(chunk,[x ' '],')'));
                        
                    end
                    
                    % convert to microns
                    x = str2double(x);
                    y = str2double(y);
                    x = x*resolution;
                    y = y*resolution;
                    
                    % ID is after coordinates, after ">", and before "</text>"
                    lbl = cell2mat(extractBetween(chunk,'>','</text>'));
                    possibleLabels{i,1} = x;
                    possibleLabels{i,2} = y;
                    possibleLabels{i,3} = lbl;
                end
            else
                disp('no labels or sites found. Site locations inferred from polygon centroids.')
            end
        end
    end
end



 %% colors
    % Prep: make a look up table of .st labels and # color codes
    % find lines that start with ".st"
    clrCdLns = find(contains(textFile{1,1},'.st'));
    
    % find the parts of those lines that correspond to style codes and color codes
    styleCodes = {NaN(numel(clrCdLns),4)};
    colorCodes = {}; colorCodeCount = 0;
    for i = 1:numel(clrCdLns)
        chunk = textFile{1,1}{clrCdLns(i),1};
        endchar = strfind(chunk,'}');
        
        % make sure we're dealing with a style definition
        if isempty(chunk) == 1
            chunk = NaN;
        elseif isempty(chunk) == 0
            
            test.none = 'fill:none;}';
            test.solid = 'fill:#';
            test.svgID = 'SVGID';
            test.swatch = 'Swatch';
            test.pattern = 'Unnamed_Pattern';
            test.stroke = 'stroke';
            test.font = 'font';
            names = fields(test);
            for j = 1:size(names,1)
                if contains(chunk, test.(names{j,:}))
                    type = names{j,:};
                    break
                end
            end
            
            switch type
                %%
                case 'solid' %extracts hex color code from the text if chunk contains it.
                    colorCodeCount = colorCodeCount+1;
                    code = cell2mat(extractBetween(chunk,'.st','{'));
                    clr = cell2mat(extractBetween(chunk,'fill:#',';'));
                    if isempty(clr)
                        code = 'NaN';
                        clr = 'NaN';
                    end
                    styleCodes{i,1} = ['st' code];
                    styleCodes{i,2} = ['#' clr];
                    styleCodes{i,3} = NaN;
                    colorCodes{colorCodeCount,1} = ['st' code];
                    colorCodes{colorCodeCount,2} = ['#' clr];
                    
                    %%
                case 'svgID'
                    code = cell2mat(extractBetween(chunk,'.st','{'));
                    clr = cell2mat(extractBetween(chunk,'url(#SVGID_','_);'));
                    if isempty(code)
                        code = 'NaN';
                    end
                    if isempty(clr)
                        clr = 'NaN';
                    end
                    styleCodes{i,1} = ['st' code];
                    styleCodes{i,3} = ['SVGID_' [clr '_']];
                    %
                    %                 code = cell2mat(extractBetween(chunk,'.st','{'));
                    %                 clr = cell2mat(extractBetween(chunk,'url(#SVGID_',');}'));
                    %                 if isempty(clr)
                    %                     code = 'NaN';
                    %                     clr = 'NaN';
                    %                 end
                    %                 styleCodes{i,1} = ['st' code];
                    %                 styleCodes{i,3} = ['SVGID_' clr];
                    
                    % look to see if SVGID is used to directly reference some styles
                    if ~isempty(find(contains(textFile{1,1},['id="SVGID_' clr '" viewBox='])))
                        % find the lines between '<pattern' and '</pattern>' nearest the use of 'id="SVGID_'
                        ptrnLn = find(contains(textFile{1,1},['id="SVGID_' clr  '_']));
                        
                        % if a string like that doesn't exist, then
                    else
                        % find where that SVGID is used
                        svgLn = find(contains(textFile{1,1},['id="SVGID_' clr '_']));
                        
                        % on the line of SVG file that defines the p
                        % find the text between '#' and '" patternTransform'
                        swtch = extractBetween(textFile{1,1}(svgLn,:),'#','" patternTransform=');
                        ptrnLn = find(contains(textFile{1,1},['id="' swtch{1,:}]));
                    end
                    
                    % find the place in the code that defines a pattern based on component style codes
                    for h = ptrnLn:size(textFile{1,1},1)
                        if contains(textFile{1,1}(h,1),'</pattern>')
                            endInd = h;
                            break
                        else
                        end
                        if h==size(textFile{1,1},1)
                            error(['the style code SVGID_' clr ' does not fit the usual pattern. You should probably get Chris to update this function to acomodate the new text patterns. Sorry :/ -Chris '])
                        end
                    end
                    bigChunk = [textFile{1,1}(ptrnLn:endInd,1)];
                    bigChunk = strjoin(bigChunk);
                    if contains(bigChunk, 'class="st')
                        % find the hex code and store it in a cell to be = the 2nd col of
                        styles = unique(extractBetween(bigChunk,'class="','" points'));
                        if isempty(styles)
                            styles = unique(extractBetween(bigChunk,'class="','" width="'));
                            if isempty(styles)
                                error('your swatches are not defined by polygons or rectangles. Modify code to accept other geometries, line 200')
                            end
                        end
                        styleCodes{i,2} = [styles]; % turn these into hashtags later
                    else
                    end
                    
                    %%
                case 'swatch'
                    code = cell2mat(extractBetween(chunk,'.st','{'));
                    clr = cell2mat(extractBetween(chunk,'url(#New_Pattern_',');}'));
                    if isempty(clr)
                        code = 'NaN';
                        clr = 'NaN';
                    end
                    styleCodes{i,1} = ['st' code];
                    styleCodes{i,2} = []; % get the hashtags
                    styleCodes{i,3} = [clr];
                    
                    ptrnLn = find(contains(textFile{1,1},['id="New_Pattern_' clr]));
                    
                    % find the place in the code that defines a pattern based on component style codes
                    for h = ptrnLn:size(textFile{1,1},1)
                        if contains(textFile{1,1}(h,1),'</pattern>')
                            endInd = h;
                            break
                        else
                        end
                        if h==size(textFile{1,1},1)
                            error(['the style code ' clr ' does not fit the usual pattern. You should probably get Chris to update this function to acomodate the new text patterns. Sorry :/ -Chris '])
                        end
                    end
                    bigChunk = [textFile{1,1}(ptrnLn:endInd,1)];
                    bigChunk = strjoin(bigChunk);
                    if contains(bigChunk, 'class="st')
                        % find the hex code and store it in a cell to be = the 2nd col of
                        styles = unique(extractBetween(bigChunk,'class="','" points'));
                        styleCodes{i,2} = [styles]; % turn these into hashtags later
                    else
                    end
                    
                    %%
                case 'pattern'
                    code = cell2mat(extractBetween(chunk,'.st','{'));
                    clr = cell2mat(extractBetween(chunk,'Unnamed_Pattern',');}'));
                    if isempty(clr)
                        code = 'NaN';
                        clr = 'NaN';
                    end
                    styleCodes{i,1} = ['st' code];
                    styleCodes{i,2} = [];
                    styleCodes{i,3} = ['Unnamed_Pattern_' clr];
                    
                    ptrnLn = find(contains(textFile{1,1},['id="Unnamed_Pattern_' clr]));
                    
                    % find the place in the code that defines a pattern based on component style codes
                    for h = ptrnLn:size(textFile{1,1},1)
                        if contains(textFile{1,1}(h,1),'</pattern>')
                            endInd = h;
                            break
                        else
                        end
                        if h==size(textFile{1,1},1)
                            error(['the style code ' clr ' does not fit the usual pattern. You should probably get Chris to update this function to acomodate the new text patterns. Sorry :/ -Chris '])
                        end
                    end
                    bigChunk = [textFile{1,1}(ptrnLn:endInd,1)];
                    bigChunk = strjoin(bigChunk);
                    if contains(bigChunk, 'class="st')
                        % find the hex code and store it in a cell to be = the 2nd col of
                        styles = unique(extractBetween(bigChunk,'class="','" points'));
                        styleCodes{i,2} = [styles]; % turn these into hashtags later
                    else
                    end
                    
                    %%
                case 'none'
                    code = cell2mat(extractBetween(chunk,'.st','{'));
                    clr = cell2mat(extractBetween(chunk,'fill:',';}'));
                    if isempty(clr)
                        code = 'NaN';
                        clr = 'NaN';
                    end
                    styleCodes{i,1} = ['st' code];
                    styleCodes{i,2} = [clr];
                    styleCodes{i,3} = NaN;
                    
                    %%
                case 'stroke'
                    code = cell2mat(extractBetween(chunk,'.st','{'));
                    clr = 'stroke';
                    if isempty(code)
                        code = 'NaN';
                        clr = 'NaN';
                    end
                    styleCodes{i,1} = ['st' code];
                    styleCodes{i,2} = [clr];
                    styleCodes{i,3} = NaN;
                    
                    %%
                case 'font'
                    code = cell2mat(extractBetween(chunk,'.st','{'));
                    clr = 'font';
                    if isempty(code)
                        code = 'NaN';
                        clr = 'NaN';
                    end
                    styleCodes{i,1} = ['st' code];
                    styleCodes{i,2} = [clr];
                    styleCodes{i,3} = NaN;
            end
        end
    end
    
    % use colorCodes to replace style codes with hashtags in styleCodes column 2
    for i = 1:size(styleCodes,1)
        if size(styleCodes{i,2},1)>1
            hashtags = {};
            hashInd = 0;
            for j = 1:size(styleCodes{i,2},1)
                clrInd = strcmp(colorCodes(:,1),styleCodes{i,2}(j));
                if sum(clrInd)==1
                    hashInd = hashInd+1;
                    hashtags{hashInd,1} = colorCodes{clrInd,2};
                elseif sum(clrInd)<1
                    % do nothing
                elseif sum(clrInd)>1
                    error('something broke in clrInd = strcmp(colorCodes(:,1),styleCodes{i,2}(j));')
                end
            end
            styleCodes{i,2} = hashtags;
        end
    end
    
    for i = 1:size(styleCodes,1)
        if iscell(styleCodes{i,2}) && size(styleCodes{i,2},1)==1
            str = cell2mat(styleCodes{i,2});
            styleCodes{i,2} = str;
        end
    end
    
%% Milad: Just in case there are no circles
    % just in case we do not have circles, we should take the coordinates
    % of the number lable (anchor point) as the center of the circle
% 
%     if isempty(circles)
%         possibleSites(:,1)=possibleLabels(:,1);
%         possibleSites(:,2)=possibleLabels(:,2);
%         possibleSites(:,3)=possibleLabels(:,3);
% 
% 
%         easyPlotSites(:,1) = possibleLabels(:,1);
%         easyPlotSites(:,2) = possibleLabels(:,2);
%     end



    %%   %% areas
    
    % find the indices of layers
    
    gStarts = find(contains(textFile{1,1},'<g id='));
    % account for nested layers (two '<g id=' right on top of each other, take away the second)
    gStarts(find(([gStarts; 0]-[0; gStarts])==1))=[];
    
    gEnds = find(contains(textFile{1,1},'</g>'));
    layers = NaN(size(gStarts,1),2);
    for i = 1:size(gStarts, 1)
        layers(i,1) = gStarts(i,:);
        try
            inds = find((gEnds>gStarts(i))&(gEnds<gStarts(i+1)));
            layers(i,2) = gEnds(inds(1));
        catch
            inds = find((gEnds>gStarts(i)));
            layers(i,2) = gEnds(inds(1));
        end
    end
    
    % find the area of S1 if it exists
    
    
    % find which layer(s) have the most polygons or paths
    for i = 1:size(layers, 1)
        bigChunks{i,1} = textFile{1,1}(layers(i,1):layers(i,2),:);
        numPolygons(i,1) = numel(find(contains(bigChunks{i,1},'<polygon class=')|contains(bigChunks{i,1},'<path class=')));
    end
    
    
    
    cellLayer = bigChunks{find(numPolygons==max(numPolygons)),1};
    
    %find all instances of "<polygon class"
    polygonLines = find(contains(cellLayer,'<polygon class=')|contains(cellLayer,'<path class='));
    
    cells = cell(size(polygonLines,1),6);
    % loop through each cell


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1:size(polygonLines,1)
        
        % find the color and area of that polygon
        chunk = cellLayer{polygonLines(i),1};
        
        % check if the line ends with '>'
        endchar = strfind(chunk,'>');
        % if not, then grab the next line down and concatenate
        loopline = 0;
        while isempty(endchar)
            loopline = loopline+1;
            chunk2 = cellLayer{polygonLines(i)+loopline,1};
            endchar = strfind(chunk2,'>');
            chunk = [chunk chunk2];
        end
        
        % find the coordinates: part of the line between points=" and "/>
        if ~(strfind(chunk, 'points')==0)
            coordsStr = cell2mat(extractBetween(chunk,'points="','"/>'));
            
            try
                if isempty(coordsStr) % check to see if the string has points
                    strings = strsplit(chunk,'points="');
                    coordsStr = strings{2};
                end
                % split up into x and y
                coordsStr = split(coordsStr);
                coordsStr(strcmp('',coordsStr)) = [];
                
                x = NaN(size(coordsStr));
                y = NaN(size(coordsStr));
                for j = 1:size(coordsStr,1)
                    str = coordsStr{j,1};
                    
                    % x is each number before the ',' and y is each number after the ','
                    coords = str2double(split(str,','));
                    
                    % convert to microns
                    coords = coords.*resolution;
                    
                    % aplit to x and y
                    x(j,1) = coords(1,1);
                    y(j,1) = coords(2,1);
                end
                
            catch % if its a path, use the helper function
                thing = extractBetween(chunk, 'd="', 'z"/>');
                thing = ['d="' cell2mat(thing) 'z"/>'];
                [xpts, ypts] = extractCoordinates([thing]);
            end
        elseif ~(strfind(chunk, 'path')==0)
            % have to use helper function
            [x, y] = extractCoordinates(chunk);
            x = x*resolution;
            y = y*resolution;
        end
        
        %find which hex codes that go with this cells style code
        styleStr = cell2mat(extractBetween(chunk,'class="','" '));
        stylCdInd = find(strcmp(styleCodes(:,1),styleStr));

        if ~isempty(stylCdInd)==1
        color = styleCodes{stylCdInd,2};
        else
         color=nan   ;
        end
        cells{i,5} = color;
        
        % package up coordinates
        vertices = {[x y]};
        cells(i,6) = vertices;
        
        
        %% if there are labels
        % find which label is inside that polygon
        % loop through that matrix of labels and test each label point to see if its in the polygon
        if exist('possibleLabels') %sum(isnan([labels{:,3}]'))==0
            try
                
                %find a site that is inside the polygon
                inSt = false(size(possibleLabels,1),1);
                for j = 1:size(possibleLabels,1)
                    xPt = possibleLabels{j,1};
                    yPt = possibleLabels{j,2};
                    labeloo = cell2mat(possibleLabels(j,3));
     
                    inSt(j) = inpolygon(xPt, yPt, [x;x(1)], [y;y(1)]); % it is a native MATLAB command
                    %
                    %                 plot(xPt, -yPt,'r*')
                    %                 plot(x, -y,'c')
                end
                
                % make a center if it doesn't exist

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 if sum(inSt)==0
%                     polyin = polyshape([x y]);
%                     [xctr,yctr] = centroid(polyin);
%                     inSt(end+1) = true; %inpolygon(xctr, yctr, x, y); !!!!!!! !!! !! wtf, this is saying its own center isn't in it
%                     %restOfSites = sites{find(in),:};
%                     % make this not overwrite everything
%                     % make this insert a cell at the right index?
%                     indEnd = size(possibleSites,1);
%                     possibleSites{indEnd+1,1} = xctr;
%                     possibleSites{indEnd+1,2} = yctr;
%                     possibleSites{indEnd+1,3} = NaN;
%                 end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % checks if in contains any matched labels to polygons.
                if sum(inSt)==0
                    close all;figure;hold on;
                    plot([x;x(1)],-[y;y(1)],'b');
                    xCtrs = [possibleLabels{:,1}]';
                    yCtrs = [possibleLabels{:,2}]';
                    plot(xCtrs,-yCtrs,'r*');
                    axis equal;
                    error('this cell is missing labels or site markers, and is ill-defined so that a center cannot be found. Re-do your SVG file!')
                elseif sum(inSt) > 1
                    error('There is more than one label or site inside a polygon. Check your SVG file!')
                end
            catch
                disp('Warning: no labels, sites, or polygon centers found. Something is seriously wrong with your SVG file!')
            end
        else
                    error('There are missing Labels' )
                    break
        end  
        

    try
        label = cell2mat(possibleLabels(inSt,3));
        cells{i,1} = label;
    catch
        cells{i,1} = 'NaNe';
    end



        %% Finding The centerpoints in the corresponding polygon
 if exist('possibleSites') %sum(isnan([labels{:,3}]'))==0
            try
                
                %find a site that is inside the polygon
                inStt = false(size(possibleSites,1),1);
                for j = 1:size(possibleSites,1)
                    xPt = possibleSites{j,1};
                    yPt = possibleSites{j,2};
                        
                    inStt(j) = inpolygon(xPt, yPt, [x;x(1)], [y;y(1)]); % it is a native MATLAB command
                    %
                    %                 plot(xPt, -yPt,'r*')
                    %                 plot(x, -y,'c')
                end
                
                % make a center if it doesn't exist
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % checks if in contains any matched labels to polygons.
                if sum(inStt)==0
                    close all;figure;hold on;
                    plot([x;x(1)],-[y;y(1)],'b');
                    xCtrs = [possibleSites{:,1}]';
                    yCtrs = [possibleSites{:,2}]';
                    plot(xCtrs,-yCtrs,'r*');
                    axis equal;
                    error('this cell is missing labels or site markers, and is ill-defined so that a center cannot be found. Re-do your SVG file!')
                elseif sum(inStt) > 1
                    error('There is more than one label or site inside a polygon. Check your SVG file!')
                end
            catch
                disp('Warning: no labels, sites, or polygon centers found. Something is seriously wrong with your SVG file!')
            end

       xctr= cell2mat(possibleSites(inStt,1));
        yctr= cell2mat(possibleSites(inStt,2));



  else
                  disp('There are missing Circles')
                    
            % if there are no labels or centers, make them from the centroids of the polygons
            polyin = polyshape([x y]);
            [xctr,yctr] = centroid(polyin);
            Centerpoint(i,1) = xctr;
            Centerpoint(i,2) = yctr;
%             Centerpoint(i,3) = 'NaNn';
%             inSt(i) = true;

 
 end

        
        
        % add site x location
        cells{i,2} = xctr;
         % add site y location
        cells{i,3} = yctr;


        % add site (LABEL)  X  location
        cells{i,10} = possibleLabels{inSt,1};
         % add site (LABEL)  y  location
        cells{i,11} = possibleLabels{inSt,2};

        
        % find the area of that polygon
        area = polyarea(x,y);
        cells{i,4} = area;
        
        % uncomment this section if you want to visually check your polygons and centers
        %             figure; hold on
        %             plot(cells{i,2}, -cells{i,3}, 'r*')
        %             plot(cells{i,6}(:,1), -cells{i,6}(:,2), 'b')
        %             close
        
    end
 %% Plot the un-chopped up map
    fig = figure;
    hold on;
    for i = 1:size(cells,1)
        plot([cells{i,6}(:,1);cells{i,6}(1,1)],-[cells{i,6}(:,2);cells{i,6}(1,2)],'color', [.7 .7 .7], 'LineWidth', 1)
        
        plot(cells{i,2},-cells{i,3},'Marker','.','color',[.7 .7 .7])  % ploting the center points
        text((cells{i,2}+1),-(cells{i,3}+1),cells{i,1},"FontSize",6)
%         plot(easyPlotSites(:,1), -easyPlotSites(:,2), 'ko');
        
    end
        if Scaling_On_Off==1

         line(Scale_info(1:2,1),Scale_info(1:2,2), "LineWidth",2, "Color",[.7 .7 .7])
         line([min(Scale_info(1:2,1))     min(Scale_info(1:2,1))],    [Scale_info(1,2)     Scale_info(3,2)] , "LineWidth",2, "Color",[.7 .7 .7])
        end

    %% Allocation of the variable
         
         Extract = cells;
   

%% Redo Identification
% We seperate the Redo sites (any numbers before or after "/"  ). We would
% use the smallest number as the the label in the  column after the
% vertices which in the current version would be  #7. Column #8 will be all
% the redo sites of the corresponding site.
% add another 2 columns as NAN or Zeros (0) to the cell/ Extract matrix 

Site_Row_Key=[0,0 ; 0,0];
kivi=1;
for J=1:1:size(Extract,1)
        Comma=0;
   k=strfind(Extract{J, 1},'/')
        if isempty (k)
           k= strfind(Extract{J, 1},',')
            Comma=1;
        end



    clear Sp1Val
   
   if ~isempty(k) % in case there are redos
            if  Comma==1
                Spl= split(Extract{J, 1},",");
            else
                
                 Spl= split(Extract{J, 1},"/");
            end

          for ji=1:1:size(Spl,1)
              Sp1Val(ji)=  str2double(Spl{ji,1});
          end

          Sorted_Redos=sort(Sp1Val,"ascend");
           Extract{J, 7}=Sorted_Redos(1);
           Extract{J, 8}=  Sorted_Redos(2:end);
            Extract{J, 9}=  Sorted_Redos(1:end);

                 
            %% This part makes a vector to map the Extract rows to each site number:
            for kik=1:1:length(Sorted_Redos)
            Site_Row_Key(kivi,1)=Sorted_Redos(kik);
            Site_Row_Key(kivi,2)=J;
            kivi=kivi+1;

            end


                                
   elseif isempty(k) % in case there are no redos
            
        Extract{J, 7}=str2double(Extract{J, 1});
        Extract{J, 8}= 0;  % there are zero redos
        Extract{J, 9}= 0;
             %% This part makes a vector to map the Extract rows to each site number:
                 Site_Row_Key(kivi,1)=str2double(Extract{J, 1});
                 Site_Row_Key(kivi,2)=J;
                 kivi=kivi+1;


   end
end

 %%  Finding the unique sites

 Unique_Sites=unique(cell2mat(Extract(:, 7))) % We find the number of unique sites and not the redos


     %% format outputs
    
    tmp=Extract;
    cells = [];
    for i = 1:size(tmp,1)
        cells(i).siteLabel = tmp{i,1};
        cells(i).siteX = tmp{i,2};
        cells(i).siteY = tmp{i,3};
        cells(i).cellArea = tmp{i,4};
        cells(i).cellColor = tmp{i,5};
        cells(i).cellVertices = tmp{i,6};
        cells(i).Siteno = tmp{i,7};
        cells(i).SiteRedo = tmp{i,8};
     end
 
     
     
     
     
     %% Patching colors
% Here is a sample to color a polygon:

% patch([cells{i,6}(:,1);cells{i,6}(1,1)],-[cells{i,6}(:,2);cells{i,6}(1,2)],'r')
% or a simpler version:

% patch([cells{i,6}(:,1)],-[cells{i,6}(:,2)],'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





  
