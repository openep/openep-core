function woi = getWindowOfInterest(userdata, varargin)
% GETWINDOWOFINTERST Returns the window of interest 
%
% Usage:
%   woi = getWindowOfInterest( userdata )
% Where:
%   userdata  - see importcarto_mem
%   woi  - two-element array specifying the window of interest relative to
%          the reference annotation
%
% GETWINDOWOFINTEREST accepts the following parameter-value pairs
%   'iEgm'    {:} | integer | array
%        - The electrogram point(s) for which the window of interst is required
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% woi = getWindowOfInterest(userdata)
% woi = getWindowOfInterest(userdata, 'iEgm', 1)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; 
iEgm = ':';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'iegm'
                iEgm = varargin{i+1};
        end
    end
end

woi = [];
if isstr(iEgm) && strcmpi(iEgm, ':')
    woi = userdata.electric.annotations.woi;
else
    for i = 1:numel(iEgm)
        woi(i,1:2) = userdata.electric.annotations.woi(iEgm(i),1:2);
    end
end

end