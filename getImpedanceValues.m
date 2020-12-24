function imp = getImpedanceValues( userdata, varargin )
% GETIMPEDANCEVALUE Returns the impedance value of given point(s)
%
% Usage:
%   imp = getImpedanceValue( userdata, varargin )
% Where:
%   userdata  - see importcarto_mem
%   imp  - the impedance values (Ohms)
%
% GETIMPEDANCEVALUES accepts the following parameter-value pairs
%   'method'    {'map'} | 'egm'
%   'points'     {':'} | int array
%   'vertices'  {':'} | int array
%
% GETIMPEDANCEVALUES Returns the impedance values. By default, impedance
% values are returned for all the points in the map. If 'method' is
% specified to be 'egm' then impedance transients are returned for each
% individual mapping point, along with time intervals for the impedances.
% If one or more 'vertices' are specified then impedance values are only 
% returned for those vertices (only valid if 'method' is 'map'). 
% If one or more 'points' is specified then impedance values are only 
% returned for those mapping points (only va;lid if 'method' is 'egm').
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% imp = getImpedanceValues(userdata, 'method', 'egm', 'points', [1 2 3]);
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 2; % UPDATE VALUE
method = 'map';
points = ':';
vertices = ':';
pointsHaveBeenSpecified = false;
verticesHaveBeenSpecified = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'points'
                points = varargin{i+1};
                pointsHaveBeenSpecified = true;
            case 'vertices'
                vertices = varargin{i+1};
                verticesHaveBeenSpecified = true;
        end
    end
end

if strcmpi(method, 'map') && pointsHaveBeenSpecified
    error('OPENEP/GETIMPEDANCEVALUES: It is not valid to specify "points" when the method is "map", use "egm" instead.');
end
if strcmpi(method, 'egm') && verticesHaveBeenSpecified
    error('OPENEP/GETIMPEDANCEVALUES: It is not valid to specify "vertices" when the method is "egm", use "map" instead.');
end

switch method
    case 'map'
        if strcmpi(vertices, ':')
            imp = userdata.surface.uni_imp_frc(:,2);
        else
            imp = userdata.surface.uni_imp_frc(vertices,2);
        end
        
    case 'egm'
        if strcmpi(points, ':')
            imp = userdata.electric.impedances;
        else
            imp.time = userdata.electric.impedances.time(points);
            imp.value = userdata.electric.impedances.value(points);
        end
end

