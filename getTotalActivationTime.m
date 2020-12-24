function tat = getTotalActivationTime( userdata, varargin )
% GETTOTALACTIVATIONTIME Returns the total activation time of the chamber
%
% Usage:
%   tat = getTotalActivationTime( userdata )
% Where:
%   userdata  - see importcarto_mem
%   tat  - the total activation time, in ms
%
% GETTOTALACTIVATIONTIME accepts the following parameter-value pairs
%   'method'    {'ptbased'} | 'ptbasedprct' | 'clinmap' | 'clinmapprct' | 'openepmap' | 'openmapprct'
%       - Specifies the method by which total activation time is calculated
%   'prct'   {2.5} | double
%       - The percentile to use for percentile mapping; only applicable if
%        'method' is one of 'ptbasedprct', 'clinmapprct' or 'openepmapprct'.
%
% GETTOTALACTIVATIONTIME Returns the total activatoin time of the chamber.
% Several alternative methods are provided, and specified by setting the
% 'method' parameter-value pair to one of the following options:
%       'ptbased'    - Calculates the difference in activation time between 
%                       the earliest and latest activation time mapping 
%                       points exported by the clinical system.
%       'ptbasedprct'- First calculates the earliest 2.5th percentile and 
%                       the latest 2.5th percentile mapping times on the 
%                       exported electrogram annotations, then calculates 
%                       the difference between the means of these sets of 
%                       activation times.
%       'clinmap'    - Calculates the difference between the earliest and 
%                       latest activation times on the local activation 
%                       time map created by the clinical mapping system
%       'clinmapprct'- First calculates the earliest 2.5th percentile and 
%                       latest 2.5th percentile mapping times on the 
%                       clinical local activation time map, then calculates 
%                       the difference between the means of these sets of 
%                       activation times.
%       'openepmap'  - Calculates the difference between the earliest and 
%                       latest activation times on the local activation 
%                       time map created by OpenEP from the exported 
%                       electrogram annotations.
%       'openepmapprct'- First calculates the earliest 2.5th percentile and 
%                       latest 2.5th percentile mapping times on the local 
%                       activation time map created by OpenEP from the 
%                       exported electrogram annotations. Then calculates 
%                       the difference between the means of these sets of 
%                       activation times.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% tat = getTotalActivationTime( userdata )
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
   error(['GETTOTALACTIVATIONTIME: Unrecognised value: ' method ' for parameter: method']); 
end
if ~isnumeric(prct)
    error(['GETEARLIESTACTIVATIONSITE: Unrecognised value ' prct ' for parameter prct.']);
end

disp(['GETTOTALACTIVATIONTIME: Using method: ' method '.']);

[~, ~, ~, tE] = getEarliestActivationSite(userdata, 'method', method, 'prct', prct);
[~, ~, ~, tL] = getLatestActivationSite(userdata, 'method', method, 'prct', prct);

tat = tL - tE;

end
