function hSurf = drawMap(userdata, varargin)
% DRAWMAP plots an OpenEP map
% Usage:
%   hSurf = drawMap(userdata)
%   hSurf = drawMap(userdata, varargin);
%
% Where:
%   hSurf - is a handle to the surface
%   userdata - is a Carto data structure
%
% DRAWMAP accepts the following parameter-value pairs
%   'data' {[]} | [d]
%       - Where d is a vector of data values and size(d) equals numel(userdata.surface.triRep.X)
%   'type'  {'act'} | 'bip' | 'force' | 'uni' | 'none' | 'cv' | 'geodesic' | 'wallthickness'
%       - Specifies type of map - activation, bipolar or unipolar voltage
%   'coloraxis' {[]} | [a b]
%       - Where a and b are real numbers. See help colorShell
%   'noLight' {false} | true
%       - If set to true no additional light will be drawn. Useful if
%       overlaying maps.
%   'usrColorMap' {[]}|cMap
%       - If set, this colormap will be used instead of the defaults
%   colorbarlocation    'north' | 'south' | 'east' | 'west' | 'northoutside' |
%                       'southoutside' | 'eastoutside' | {'westoutside'}
%   'orientation' {'AP'} | 'PA'
%       - Specifies the view as AP or PA. LAO, RAO, LL, RL yet to be
%       defined
%   'colorfillthreshold'   {10} | c
%       - Where c is a scalar value; defaulting to 10mm
%
% DRAWMAP is a wrapper function for colorShell.m which allows OpenEP data 
% to be plotted. 
%
% Author: Steven Williams (2016) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------
% hSurf = drawMap(userdata, 'type', 'act')
% set(hSurf, 'facealpha', .5);
% hSphere = plot3dsphere(userdata.electric.egmSurfX, 'r', 3, 16, 'stepthrough', true)
% ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% Parse the input variables
p = inputParser;
p.addRequired('userdata');
p.addParameter('data', [], @isnumeric);
p.addParameter('type', 'act', @ischar);
p.addParameter('coloraxis', [], @isnumeric);
p.addParameter('noLight', false, @isbool);
p.addParameter('usrColorMap', [], @isnumeric);
p.addParameter('colorbarlocation', 'westoutside', @(x)ischar(x) && (strcmpi(x, 'north') ...
    || strcmpi(x, 'south') || strcmpi(x, 'east') || strcmpi(x, 'west') ...
    || strcmpi(x, 'northoutside') || strcmpi(x, 'southoutside') ...
    || strcmpi(x, 'eastoutside') || strcmpi(x, 'westoutside') ...
    ));
p.addParameter('orientation', 'ap', @ischar);
p.addParameter('colorfillthreshold', 10, @isnumeric);
p.parse(userdata, varargin{:});
inputs = p.Results;
userdata = inputs.userdata;
DATAMANUAL = inputs.data;
type = inputs.type;
coloraxis = inputs.coloraxis;
noLight = inputs.noLight;
usrColorMap = inputs.usrColorMap;
colorbarlocation = inputs.colorbarlocation;
orientation = inputs.orientation;
DISTANCETHRESH = inputs.colorfillthreshold;

if ~any(strcmpi(type, {'act', 'bip', 'force', 'uni', 'none', 'cv', 'geodesic', 'wallthickness'}))
    error('DRAWMAP: Invalid parameter-value pair; type must be one of act, bip, force, uni, none, cv or geodesic')
end


% Draw the surface
hSurf = trisurf(getMesh(userdata), 'edgecolor', 'none');
axis equal vis3d
set(hSurf, 'facecolor', [.5 .5 .5]);

% Draw the free boundary, i.e. valve
drawFreeBoundary(getMesh(userdata), [0 0 0]);

% Set up variables for colorShell
switch type
    case 'act'
        DATATYPE = 'activation';
        if ~isempty(DATAMANUAL)
            DATA = DATAMANUAL;
        else
            DATA = userdata.surface.act_bip(:,1);
        end
    case 'bip'
        DATATYPE = 'voltage';
        if ~isempty(DATAMANUAL)
            DATA = DATAMANUAL;
        else
            DATA = userdata.surface.act_bip(:,2);
        end
    case 'force'
        DATATYPE = 'force';
        if ~isempty(DATAMANUAL)
            DATA = DATAMANUAL;
        else
            DATA = userdata.surface.uni_imp_frc(:,3);
        end
    case 'uni'
        DATATYPE = 'voltage';
        if ~isempty(DATAMANUAL)
            DATA = DATAMANUAL;
        else
            DATA = userdata.surface.uni_imp_frc(:,1);
        end
    case 'cv'
        DATATYPE = 'cv';
        if ~isempty(DATAMANUAL)
            DATA = DATAMANUAL;
        else
            DATA = getConductionVelocity(userdata);
        end
    case 'geodesic'
        DATATYPE = 'geodesic';
        if ~isempty(DATAMANUAL)
            DATA = DATAMANUAL;
        else
            % by default draw the geodesics from the first electrode positon
            distances = getData2MeshGeodesicDistances(userdata);
            DATA = NaN(size(getVertices(userdata),1),10);
            int = openEpDataInterpolator();
            for i = 1:20%size(distances,1)
                DATA(:,i) = int.interpolate(getVertices(userdata,'used',true), distances(i,:)', getVertices(userdata,'used',false));
            end
        end
    case 'wallthickness'
        DATATYPE = 'wallthickness';
        if ~isempty(DATAMANUAL)
            DATA = DATAMANUAL;
        else
            error('OPENEP/drawMap: automatic calculation of wall thickness data not supported')
        end
end

if ~strcmpi(type, 'none')
    if isempty(coloraxis)
        t_min = min(DATA);
        t_max = max(DATA);
    else
        t_min = coloraxis(1);
        t_max = coloraxis(2);
    end
    
    % Call colorShell
    if ~all(isnan(DATA))
        colorShell(hSurf, userdata.electric.egmSurfX, DATA, DISTANCETHRESH ...
            , 'showcolorbar', 'show' ...
            , 'coloraxis', [t_min t_max] ...
            , 'datatype', DATATYPE ...
            , 'interpolation', 'off' ...
            , 'colorbarlocation', 'north' ...
            , 'usrColorMap', usrColorMap ...
            , 'colorbarlocation', colorbarlocation ...
            );
    else
        disp('no data to color the shell with')
    end    
end

% Adjust the light/material
material dull
if ~noLight
    cameraLight;
end

% Adjust the viewpoint
switch lower(orientation)
    case 'ap'
        set(gca, 'cameraposition', get(gca, 'cameratarget') + [0 0 700])
    case 'pa'
        set(gca, 'cameraposition', get(gca, 'cameratarget') - [0 0 700])
end
set(gca, 'cameraupvector', [0 1 0]);

% Adjust the appearances
set(gcf, 'color', 'white');
axis off;


