function interpData = generateInterpData(userdata, datatype, varargin)
% GENERATEINTERPDATA Performs spatial interpolation of scalar data
%
% Usage:
%   interpData = generateInterpData(userdata, datatype)
% Where:
%   userdata - see, importcarto_mem
%   datatype - the desired data type to return
%   interpData - is the interpolated data
%
% GENERATEINTERPDATA accepts the following parameter-value pairs
%   'interMethod'    'nearest'|'linear'|{'natural'}
%       - The interpolation method, default to natural
%   'exterMethod'    {nearest}|linear|none
%       - The extrapolation method, defaults to linear
%   'distanceThresh' {10}|double
%       - The distance threshold, d, default 10mm
%
% GENERATEINTERPDDATA performs spatial interpolation of scalar data.
% Userdata and datatype are mandatory arguments. Datatype may be one of:
%         'bip-map' - bipolar voltage; from the exported voltage values
%         'uni-map' - unipolar voltage; from the exported voltage values
%         'lat-map' - local activation time; from the annotated electrograms
%         'bip-egm' - bipolar voltage; measured by OpenEP on the egms (NOT IMPLEMENTED)
%         'uni-egm' - unipolar voltage; measured by OpenEP on the egms (NOT IMPLEMENTED)
%         'lat-egm' - local activation time; measured by OpenEP on the egms (NOT IMPLEMENTED)
%         'cv' - conduction velocity
%         'egmduration' - electrogram duration
% GENERATEINTERPDATA removes any NaN values in data (and their
% corresponding location(s) in coords) before calling scatteredInterpolant
% with the interpolation/extrapolation methods specified. Any values greater
% than distancethresh are removed.
%
% TODO: Separate the interpolation function into a subroutine to permit extensibility
%
% Author: Steven Williams (2018) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% interpData = generateInterpData(userdata, 'bip-map');
% interpData = generateInterpData(userdata, 'lat-map');
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------
     
nStandardArgs = 2; 
interMethod = 'natural';
exterMethod = 'nearest';
distanceThresh = 10;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'intermethod'
                interMethod = varargin{i+1};
                disp(interMethod);
            case 'extermethod'
                exterMethod = varargin{i+1};
                disp(exterMethod);
            case 'distancethresh'
                distanceThresh = varargin{i+1};
        end
    end
end
if ~any(strcmpi(interMethod, {'nearest' 'linear' 'natural'}))
    error(['GENERATEINTERPDATA: Unrecognised value ' interMethod ' for parameter interMethod']);
end
if ~any(strcmpi(exterMethod, {'nearest' 'linear' 'none'}))
    error(['GENERATEINTERPDATA: Unrecognised value ' exterMethod ' for parameter exterMethod']);
end
if ~any(strcmpi(datatype, {'bip-map' 'uni-map' 'lat-map' 'bip-egm' 'uni_egm' 'lat-egm' 'cv' 'egmduration'}))
    error(['GENERATEINTERPDATA: Unrecognised value ' datatype ' for parameter datatype']);
end

% Get the vertices of the trirep
pts = getMesh(userdata).X;

% Get the data
switch lower(datatype)
    case 'bip-map'
        data = userdata.electric.voltages.bipolar;
        coords = userdata.electric.egmX;
    case 'uni-map'
        data = userdata.electric.voltages.unipolar;
        coords = userdata.electric.egmX;
    case 'bip-egm'
        warning('TODO: bip-egm to be implemented')
        return;
        % to be implemented
    case 'uni-egm'
        warning('TODO: uni-egm to be implemented')
        return;
        % to be implemented
    case 'cv'
        warning('TODO: cv to be implemented')
        return;
        % to be implemented
    case 'lat-map'
        data = userdata.electric.annotations.mapAnnot - userdata.electric.annotations.referenceAnnot;
        coords = userdata.electric.egmX;
    case 'egmduration'
        data = getElectrogramDuration(userdata);
        coords = userdata.electric.egmX;
end

% Remove any points that are annotated outwith the region of interest
iVp = getMappingPointsWithinWoI(userdata);
data(~iVp) = NaN;
coords(~iVp,:) = NaN;

% Do the interpolation, once for each column in data
interpData = NaN(length(pts),size(data,2));
for i = 1:size(data,2)
    % Remove any data with NaNs
    tempData = data(:,i);
    tempCoords = coords;
    iNaN = isnan(tempData);
    tempData(iNaN) = [];
    tempCoords(iNaN,:) = [];
    
    % Interpolate the data
    warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId');
    if numel(tempData)>=4
        % with more than 3 data points we can do trianulation-based interpolation
        F = scatteredInterpolant(tempCoords(:,1), tempCoords(:,2), tempCoords(:,3) ...
            , tempData ...
            , interMethod ... % interpolation {nearest|linear|natural}
            , exterMethod ... % extrapolation {nearest|linear|none}
            );
        interpData(:,i) = F(pts);
        
    elseif numel(tempData)>=3
        % do manual nearest neighbour interpolation
        tempData2 = data(:,i);
        warning('GENERATESURFACECDATA: Less than 4 data points, performing nearest neighbour interpolation');
        [vertices, ~] = findclosestvertex(coords, pts);
        interpData(:,i) = tempData2(vertices);
    else
        interpData(:,i) = NaN;
    end
    
    % work out if there are any points on the surface that are less than 0
    interpData(interpData(:,i)<0,i) = 0;
    
    % work out which points on the surface are too far away from real data
    if ~isempty(tempCoords)
        id = knnsearch(tempCoords, pts);
        cPts = tempCoords(id,:); %c for closest
        d = distBetweenPoints(cPts, pts);
        thresholdDistance = zeros(size(d));
        thresholdDistance(d>distanceThresh) = 1;
        nanset = logical(thresholdDistance);
        interpData(nanset,i) = NaN;
    end
    
    % Remove any data which is associated with vertices that are not 
    % referenced by the triangulation
    [~,isVertUsed] = repack(getMesh(userdata, 'type', 'triangulation'));
    interpData(~isVertUsed,i) = NaN;
end

end

