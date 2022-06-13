function [dataField, hColBar] = colorShell(hSurf, pts, data, t, varargin)
% COLORSHELL shades the surface hSurf with data at pts
%
% Usage:
%   colorShell(hSurf, pts, data, t)
%   colorShell(..., 'propertyname', propertyvalue)
%
% Where:
%   hSurf - see plotVelocityGeometry
%   pts - the data position coordinates. Size mpts x ndim
%   data - the data at each location of pts. Size mpts x ndata
%   t - threshold distance from pts greater than which shell will not be 
%       coloured
%
%   dataField - is the data which has been drawn
%   hColBar - is a handle to the color bar
%
% COLORSHELL shades the surface of a shell according to the scalar field
% data at location pts.  If the property-value pair 'coloraxis' is set then 
% color thresholding is performed:
%   triangles with |data|>max(caxis) - purple
%   triangles with |data|<min(caxis) - brown
% The property-value pair datatype is used to tell the colormap whether to
% flip or not (otherwise jet is back-to-front for isochronal maps)
%
% Properties:
%   showcolorbar        'show' | {'hide'}
%   colorbarlocation    'north' | 'south' | 'east' | {'west'} | 'northoutside' |
%                       'southoutside' | 'eastoutside' | 'westoutside'
%   coloraxis           [caxismin caxismax] | {[min(min(data)) max(max(data))]}
%   datatype            {'activation'} | 'force' | 'voltage'
%   interpolation       {'on'} | 'off'
%   nanset              logical array indexing into hSurf points
%   faceAlpha           integer
%   usrColorMap         color map matrix
%
% Author: Steven Williams (2013)
% Modifications -
%   Improved color thresholding - Steven Williams (2016)
%
% Info on Code Testing:
% ---------------------------------------------------------------
% % 1. Color by height
%     DISTANCETHRESH = 12;
%     [elec, hSurf, ~, ~] = plotVelocityGeometry(VelocityData, 'all', 'bip');
%     data = elec.Position(:,3);
%     colorShell(hSurf, elec.Position, data, DISTANCETHRESH);
%
% % 2. Color by a random dataset
%     DISTANCETHRESH = 12;
%     [elec, hSurf, ~, ~] = plotVelocityGeometry(VelocityData, 'all', 'bip');
%     data = rand(size(elec.Position(:,3)));
%     colorShell(hSurf, elec.Position, data, DISTANCETHRESH);
%
% % 3. Color by a single site with high values
%     SITE = 4;
%     DISTANCETHRESH = 12;
%     [elec, hSurf, ~, ~] = plotVelocityGeometry(VelocityData, 'all', 'bip');
%     data = zeros(size(elec.Position(:,3)));
%     data(1+10*SITE:10+10*SITE) = 1;
%     colorShell(hSurf, elec.Position, data, DISTANCETHRESH);
%
% % 4. Color by a number of different sites with high values and step thru
%     SITE = [1 2 3 4 5 6 7 8 9 10];
%     DISTANCETHRESH = 12;
%     [elec, hSurf, ~, ~] = plotVelocityGeometry(VelocityData, 'all', 'bip');
%     data = zeros(size(elec.Position(:,3)));
%     data = repmat(data, 1, numel(SITE));
%     for i = 1:numel(SITE)
%         data(10*(SITE(i)-1)+1:10+10*(SITE(i)-1), i) = 1;
%     end
%     colorShell(hSurf, elec.Position, data, DISTANCETHRESH);
%
% % 5. Color by individual points with high values and step thru
%     SITE = 1:100;
%     DISTANCETHRESH = 12;
%     [elec, hSurf, ~, ~] = plotVelocityGeometry(VelocityData, 'all', 'bip');
%     data = zeros(size(elec.Position(:,3)));
%     data = repmat(data, 1, numel(SITE));
%     for i = 1:numel(SITE)
%         data(SITE(i), i) = 1;
%     end
%     colorShell(hSurf, elec.Position, data, DISTANCETHRESH);
%
% % 6. Force shell
%     clf;
%     hSurf = trisurf(userdata.surface.triRep);
%     axis equal vis3d;
%     DISTANCETHRESH = 4;
%     colorShell(hSurf, userdata.electric.egmSurfX, userdata.electric.force.force, DISTANCETHRESH ... 
%       , 'showcolorbar', 'show' ...
%       , 'coloraxis', [1 40] ...
%       , 'datatype', 'force' ...
%     );
%
% % 7. Activation map based on interpolated Carto data
%     clf
%     hSurf = trisurf(userdata.surface.triRep);
%     axis equal vis3d
%     DISTANCETHRESH = 10;
%     t_min = min(userdata.surface.act_bip(:,1));
%     t_max = max(userdata.surface.act_bip(:,1));
%     colorShell(hSurf, userdata.electric.egmSurfX, userdata.surface.act_bip(:,1), DISTANCETHRESH ...
%         , 'showcolorbar', 'show' ...
%         , 'coloraxis', [t_min t_max] ...
%         , 'datatype', 'activation' ...
%         , 'interpolation', 'off' ...
%         );
%
% % 8. Color by face data values from VKT file
%     tr = hVtk.getTriRep
%     hSurf = trisurf(tr)
%     [~, allCentroids] = tricentroid(tr);
%     faceData = hVtk.CellData{1};
%     colorShell(hSurf, allCentroids, faceData, Inf);
%     axis equal vis3d
%     hVtk.CellDataNames{1}
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% Validate input
p = inputParser;
p.addRequired('hSurf', @ishandle);
p.addRequired('pts', @isnumeric);
p.addRequired('data', @isnumeric);
p.addRequired('t', @isnumeric);
p.addParameter('showcolorbar', 'hide', @(x)ischar(x) &&(strcmpi(x,'hide') || strcmpi(x, 'show')));
p.addParameter('colorbarlocation', 'west', @(x)ischar(x) && (strcmpi(x, 'north') ...
    || strcmpi(x, 'south') || strcmpi(x, 'east') || strcmpi(x, 'west') ...
    || strcmpi(x, 'northoutside') || strcmpi(x, 'southoutside') ...
    || strcmpi(x, 'eastoutside') || strcmpi(x, 'westoutside') ...
    ));
p.addParameter('coloraxis', [], @isnumeric);
p.addParameter('nanset', [], @islogical);
p.addParameter('datatype', 'activation', @isstr);
p.addParameter('colororder', 'forward', @isstr); 
p.addParameter('interpolation', 'on', @isstr);
p.addParameter('faceAlpha', 1, @isnumeric);
p.addParameter('usrColorMap', [], @isnumeric);
p.parse(hSurf, pts, data, t, varargin{:});
inputs = p.Results;
showcolorbar = inputs.showcolorbar;
colorbarlocation = inputs.colorbarlocation;
coloraxis = inputs.coloraxis;
datatype = inputs.datatype;
doInterpolation = inputs.interpolation;
nanset = inputs.nanset;
faceAlpha = inputs.faceAlpha;
usrColorMap = inputs.usrColorMap;

% Increase pts to three dimensions if necessary
if size(pts,2) == 2
    pts(:,3) = 0;
end

% Get a handle to the figure and axis
hAx = get(hSurf, 'parent');
hFig = getParentFigure(hSurf);

% Draw the sliderbar, if there is more than one dataset
if size(data,2)>1
    hSl = uicontrol(hFig, 'style', 'slider' ...
        , 'units', 'normalized' ...
        , 'position', [.8 .1 .15 .05] ...
        , 'callback', @adjustSlider ...
        , 'min', 1 ...
        , 'max', size(data,2) ...
        , 'value', 1 ...
        , 'sliderstep', [1/(size(data,2)-1) 1/(size(data,2)-1)] ...
        );
end

% Call adjust slider for the first time
adjustSlider();

    function adjustSlider(~, ~, ~)
        % Get the value of the sliderbar, if it exists
        if exist('hSl', 'var')
            i = round(get(hSl, 'value'));
        else
            i = 1;
        end
        
        % Every point in the surface - r for required.
        rPts = get(hSurf, 'Vertices');
        
        % Perform interpolation
        switch doInterpolation
            case 'on'
                % Set up the interpolation
                warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
                tempdata = data(:,i);
                temppts = pts;
                indexRemove = isnan(tempdata);
                tempdata(indexRemove,:)=[];
                temppts(indexRemove,:)=[];
                F = scatteredInterpolant(temppts, tempdata ...
                    , 'natural' ... %interpolation method - linear, nearest or natural
                    , 'linear' ... %extrapolation method - linear, nearest or none
                    ); 
                dataField = F(rPts);
            case 'off'
                dataField = data(:,i);
                temppts = pts;
        end
        
        if isempty(dataField)
            return
        end
        
        % Find the shortest distance from every point in the surface to the dataset
        % pointcloud
        id = knnsearch(temppts, rPts);
        cPts = temppts(id,:); %c for closest
        d = distBetweenPoints(cPts, rPts);
        
        % if coloraxis has been specified set up for raw data thresholding now
        tAboveCAxis = zeros(size(d));
        tBelowCAxis = zeros(size(d));
        if ~isempty(coloraxis)
            tAboveCAxis(dataField>max(coloraxis)) = 1;
            tBelowCAxis(dataField<min(coloraxis)) = 1;
        end
        
        % Define magenta and brown
        magenta = [1 0 1] * .7;
        red = [0.8 0 0];
            
        % Get the appropriate color map
        switch datatype
            case 'activation'
                colormap(ones([64,3]));
                cMap = hsv;
                cMap(end-8:end,:) = [];
                cBarTitle = 'LAT (ms)';
                magenta = cMap(end,:);
                red = cMap(1,:);
                
            case 'force'
                cMap = jet;
                cBarTitle = 'Force (g)';
                
            case 'voltage'
                cMap = flipud(jet);
                cMap(1:length(cMap)/10,:) = [];
                cBarTitle = 'Voltage (mV)';
                magenta = [1 0 1] * .8;
                red = cMap(1,:);
                %red = [1 0 0];
                
            case 'cv'
                cMap = colormap(hot);
                cBarTitle = 'CV (m/s)';
                magenta = cMap(end,:);
                red = cMap(1,:);
                
            case 'geodesic'
                cMap = flipud(colormap(winter));
                cBarTitle = 'Distance (cm)';
                
            case 'wallthickness'
                cMap = flipud(jet);
                cMap(1:length(cMap)/10,:) = [];
                magenta = [1 0 1] * .8;
                red = cMap(1,:);
                cBarTitle = 'Wall Thickness (mm)';

            case 'weights'
                cMap = colormap(parula);
                cBarTitle = 'Normalised distance weights';

            case 'labels'
                cBarTitle = 'labels';

        end
        if ~isempty(usrColorMap)
            cMap = usrColorMap;
        end

        % Limit the size of the colormap
        if length(cMap) > 256
            cMap = downsample(cMap, round(length(cMap)/256));
        end
        %colormap(hAx, cMap);
        disp(['Length of cMap is: ' num2str(length(cMap))]);
        
        scaledDataField = round(scaleData( dataField(~tAboveCAxis & ~tBelowCAxis), [1 size(cMap,1)] ));
        iC = ones(size(dataField));
        iC(find(~tAboveCAxis & ~tBelowCAxis)) = scaledDataField; %#ok<FNDSB>
        thresholdNaN = isnan(iC);
        iC(thresholdNaN) = 1;
        col = cMap(iC,:);
        
        % At this point the color map col is correct for all non-nan points
        % p satisfying the expression min(caxis)<dataField(p)<max(caxis)
        
        % Threshold based on coloraxis
        if ~isempty(coloraxis)
            col(logical(tAboveCAxis),:) = repmat(magenta, [numel(find(tAboveCAxis)),1]);
            col(logical(tBelowCAxis),:) = repmat(red, [numel(find(tBelowCAxis)),1]);
        end
        
        % Threshold based on the distance t
        thresholdDistance = zeros(size(d));
        thresholdDistance(d>t) = 1;
        
        % Create the color data matrix, making any vertex grey when it is greater
        % than t distance from the point cloud.
        col(logical(thresholdDistance),:) = .6; %ie [.6 .6 .6] grey
        col(logical(thresholdNaN),:) = .6; % ie [1 1 1] white
        col(logical(nanset),:) = .6;
        
        % Update the color of the surface
        set(hSurf, 'FaceVertexCData', col ...
            ,'FaceColor', 'interp' ...
            ,'FaceAlpha', faceAlpha ...
            );
        
        % Set the figure color axis
        % colorAxis represents the range of the whole data
        % coloraxis represents the range specified by the user
        colorAxis = [min(min(data)) max(max(data))];
        caxis(colorAxis);
        if ~isnan(colorAxis(1)) && ~isnan(colorAxis(2))
            % Create the figure color map
            if range(coloraxis)>0
                % a range of thresholds is specified
                incrementOfCMap = range(coloraxis)/(length(cMap)-1);
                nMagenta = round((colorAxis(2) - coloraxis(2))/incrementOfCMap);
                nBrown = round((coloraxis(1) - colorAxis(1))/incrementOfCMap);
                cMapMagenta = repmat(magenta, [nMagenta, 1]);
                cMapRed = repmat(red, [nBrown, 1]);
                figureCMap = vertcat(cMapRed, cMap, cMapMagenta);
            else
                % there is a single threshold
                thresh = min(coloraxis);
                nColors = 256;
                nMagenta = ceil( ((max(data(:,i)) - thresh)/range(data(:,i))) * nColors );
                nBrown = ceil( ((thresh)/range(data(:,i))) * nColors );
                cMapMagenta = repmat(magenta, [nMagenta, 1]);
                cMapRed = repmat(red, [nBrown, 1]);
                figureCMap = vertcat(cMapRed, cMapMagenta);
            end
            colormap(figureCMap);
        end
        
        % Show the colorbar
        if strcmpi(showcolorbar, 'show')
            % show the color bar
            colorbar('peer', hAx, 'delete');
            hColBar = colorbar('peer', hAx, 'location', colorbarlocation);
            hColBar.Label.String = cBarTitle;
        elseif strcmpi(showcolorbar, 'hide')
            % hide the color bar
            colorbar('peer', hAx, 'delete');
        end
        
        % Return the scaledData if asked for
        if strcmpi(datatype,'force')
            dataField(logical(thresholdDistance),:) = 0;
            dataField(logical(thresholdNaN),:) = 0;
        end
        if strcmpi(datatype,'activation')
            dataField(logical(thresholdDistance),:) = NaN;
        end
        
    end %adjustSlider()
end %colorShell()