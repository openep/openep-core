function [X, surfX, iPoint, t] = getEarliestActivationSite( userdata, varargin )
% GETEARLIESTACTIVATIONSITE Returns the earliest activation site
%
% Usage:
%   [X] = getEarliestActivationSite( userdata )
%   [X, surfX] = getEarliestActivationSite( userdata )
%   [X, surfX, iPoint] = getEarliestActivationSite( userdata )
%   [X, surfX, iPoint, t] = getEarliestActivationSite( userdata )
%
% Where:
%   userdata    - see importcarto_mem
%   X           - the earliest activation site, in x,y,z format. For
%                   map-based methods (i.e. `clinmap`, `clinmapprct`,
%                   `openepmap` and `openepmapprct`), X is identical to
%                   surfX.
%   surfX       - the surface projection of the eariest activation site
%   iPoint      - the closest mapping point to the earliest activation
%                   site. For point-based methods (i.e. `ptbased` or 
%                   `ptbasedprct`), iPoint indexes into userdata.electric.
%                   For map-based methods (i.e. `clinmap`, `clinmapprct`,
%                   `openepmap`, `openepmapprct`), iPoint indexes into
%                   getMesh(userdata).X. For percentile methods (i.e.
%                   `ptbasedprct`, `cinmapprct` or `openepmapprct`) iPoint
%                   returns all the points that were identified within the
%                   relevant percentile.
%   t           - the calculated earliest activation time, relative to the
%                   reference annotation
%
% GETEARLIESTACTIVATIONSITE accepts the following parameter-value pairs
%   'method'    {'ptbased'} | 'ptbasedprct' | 'clinmap' | 'clinmapprct' | 'openepmap' | 'openmapprct'
%       - Specifies the method by which the earliest activation is
%           calculated.
%   'prct'       {2.5} | double
%       - The percentile to use for percentile mapping; only applicable if
%       `'method'` is `'ptbasedprct'`, `'clinmapprct'` or `'openmapprct'`
%
% GETEARLIESTACTIVATIONSITE By identifying the earliest activatin site,
% this function can be used, for example, to identify the site in the
% chamber closest to the pacing site. Several alternative methods are provided for caluclating the latest 
% activation site, specified by setting the 'method' parameter-value pair 
% to one of the following options:
%       `ptbased`    - Calculates the earliest activation time using 
%                       mapping points exported by the clinical system.
%       `ptbasedprct`- Calculates the 0-2.5th percentile mapping 
%                       times on the exported electrogram annotations, then 
%                       calculates the mean of this sets of activation times.
%       `clinmap`    - Calculates the earliest activation time on the 
%                       local activation time map created by the clinical 
%                       mapping system
%       `clinmapprct`- First calculates the 0-2.5th percentile 
%                       mapping times on the clinical local activation time 
%                       map, then calculates the mean of these sets of 
%                       activation times.
%       `openepmap`  - Calculates the earliest activation time on the local 
%                       activation time map created by OpenEP from the 
%                       exported electrogram annotations.
%       `openepmapprct`- First calculates the 0-2.5th percentile 
%                       mapping times on the local activation time map 
%                       created by OpenEP from the exported electrogram 
%                       annotations. Then calculates the mean of this set of 
%                       activation times.
%
% See also GETTOTALACTIVATIONTIME, GETLATESTACTIVATIONSITE.
%
% Author: Steven Williams (2019) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%   Steven Williams 2020: Added the additional methods
%
% Info on Code Testing:
% ---------------------------------------------------------------
% [X, surfX, iPoint, t] = getEarliestActivationSite( userdata );
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
method = 'ptbased';
prct = 2.5;
if nargin > nStandardArgs
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'prct'
                prct = varargin{i+1};
        end
    end
end
if ~any(strcmpi(method, {'ptbased' 'ptbasedprct' 'clinmap' 'clinmapprct' 'openepmap' 'openepmapprct'}))
   error(['GETEARLIESTACTIVATIONSITE: Unrecognised value: ' method ' for parameter: method']); 
end
if ~isnumeric(prct)
    error(['GETEARLIESTACTIVATIONSITE: Unrecognised value ' prct ' for parameter prct.']);
end

%new version
iP = getMappingPointsWithinWoI(userdata);
mapAnnot = userdata.electric.annotations.mapAnnot;
mapAnnot(~iP) = NaN;
referenceAnnot = userdata.electric.annotations.referenceAnnot;

switch method
    case 'ptbased'
        % find the earliest point on the map, defined as the earliest map
        % annotation within the window of interest
        [minMapAnnot,iPoint] = nanmin(mapAnnot);
        X = userdata.electric.egmX(iPoint,:);
        surfX = userdata.electric.egmSurfX(iPoint,:);
        t = minMapAnnot - referenceAnnot(iPoint,1);
    case 'ptbasedprct'
        % find the earliest points on the map and average their locations
        [~, I] = sort(mapAnnot, 'ascend');
        iPoint = I(1:floor(numel(mapAnnot)*prct/100));
        X = mean(userdata.electric.egmX(iPoint,:));
        surfXtemp = mean(userdata.electric.egmSurfX(iPoint,:));
        [surfXtr, distances] = findclosestvertex(getMesh(userdata), surfXtemp);
        surfX = getMesh(userdata).X(surfXtr,:);
        t = mean(mapAnnot(iPoint) - referenceAnnot(iPoint));
        disp(['Nearest surface point to [' num2str(surfXtemp) '] is [' num2str(surfX) '], distance: ' num2str(distances) '.'])
    case 'clinmap'
        % 1. Get clinical activation data
        [t,iPoint] = nanmin(userdata.surface.act_bip(:,1));
        X = getMesh(userdata).X(iPoint,:);
        surfX = X;
    case 'clinmapprct'
        act = userdata.surface.act_bip(:,1);
        [~, I] = sort(act, 'ascend');
        iPoint = I(1:floor(numel(act)*prct/100));
        X = mean(getMesh(userdata).X(iPoint,:));
        surfXtr = findclosestvertex(getMesh(userdata),X);
        surfX = getMesh(userdata).X(surfXtr,:);
        t = mean(act(iPoint));
    case 'openepmap'
        % 1. Get OpenEP interpolated activation data
        act = generateInterpData(userdata, 'lat-map');
        [t,iPoint] = nanmin(act);
        X = getMesh(userdata).X(iPoint,:);
        surfX = X;
    case 'openepmapprct'
        % 1. Get OpenEP interpolated activation data
        act = generateInterpData(userdata, 'lat-map');
        [~, I] = sort(act, 'ascend');
        iPoint = I(1:floor(numel(act)*prct/100));
        X = mean(getMesh(userdata).X(iPoint,:));
        surfXtr = findclosestvertex(getMesh(userdata),X);
        surfX = getMesh(userdata).X(surfXtr,:);
        t = mean(act(iPoint));
end


end
