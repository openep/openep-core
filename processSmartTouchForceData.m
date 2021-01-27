function hSurf = processSmartTouchForceData(userdata, varargin)
% PROCESSSMARTTOUCHFORCEDATA Gets a shell with force data
% Usage:
%   [tr f] = processSmartTouchForceData(userdata)
% Where:
%   userdata - is a Carto dataset
%   f - are the 1000ms window forces
%   loc - are the coordinates of the data in f (size: length(f) * 3) 
%
% PROCESSSMARTTOUCHFORCEDATA parses the force data at ablation points from
% userdata.
%
% Properties:
%   drawforcemap:       {'false'} | 'true'
%   * distancethreshold:  numeric, default=4
%   * showcolorbar            show | {hide}
%   * colorbarlocation    north | south | east | {west} | northoutside |
%                       southoutside | eastoutside | westoutside
%   * coloraxis           [caxismin caxismax] | {[min(min(data)) max(max(data))]}
% * Passed directly to colorShell, and ignored if drawforcemap='false' or
% empty
%
% Author: Steven Williams (2013) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% processSmartTouchForceData(userdata, 'showcolorbar', 'show', 'coloraxis', [10 20])
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% Validate input
p = inputParser;
p.addRequired('userdata', @(x)isfield(x, 'cartoFolder'));
p.addParamValue('distancethreshold', 4, @isnumeric);
p.addParamValue('showcolorbar', 'hide', @(x)isstr(x) &&(strcmpi(x,'hide') || strcmpi(x, 'show'))); %#ok<DISSTR>
p.addParamValue('colorbarlocation', 'west', @(x)isstr(x) && (strcmpi(x, 'north') ...
    || strcmpi(x, 'south') || strcmpi(x, 'east') || strcmpi(x, 'west') ...
    || strcmpi(x, 'northoutside') || strcmpi(x, 'southoutside') ...
    || strcmpi(x, 'eastoutside') || strcmpi(x, 'westoutside') ...
    )); %#ok<DISSTR>
p.addParamValue('coloraxis', [], @isnumeric);
p.parse(userdata, varargin{:});
inputs = p.Results;
userdata = inputs.userdata;
distanceThreshold = inputs.distancethreshold;
showcolorbar = inputs.showcolorbar;
colorbarlocation = inputs.colorbarlocation;
coloraxis = inputs.coloraxis;

% Get the ablation data
f = userdata.electric.force.force;
f(~strcmpi(userdata.electric.tags,'ablation'))=[];

% Get the ablation coordinates
loc = userdata.electric.egmSurfX(strcmpi(userdata.electric.tags,'ablation'),:);
if isempty(loc)
    error('PROCESSSMARTTOUCHFORCEDATA: No ablation points found; are you sure you have loaded the correct map?');
end

% Draw the force map
hSurf = trisurf(userdata.surface.triRep);
axis equal vis3d;
scaledData = colorShell(hSurf, loc, f, distanceThreshold ...
        , 'datatype', 'force' ...
        , 'coloraxis', coloraxis ...
        , 'showcolorbar', showcolorbar ...
        , 'colorbarlocation', colorbarlocation);

% Ask the user if they want to save the shell as a VTK file
[filename, pathname, filterindex] = uiputfile('*.vtk', 'Save the shell as a VTK file?');
if filterindex~=0
    % The user wants to save the file
    writeTriRep2VTK(userdata.surface.triRep, scaledData, 'outputfile', [pathname filesep() filename]);
end
end
