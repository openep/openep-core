function map = getMapAnnotation(userdata, varargin)
% GETMAPANNOTATION Returns the value(s) of the map annotation 
%
% Usage:
%   map = getMapAnnotation( userdata )
% Where:
%   userdata  - see importcarto_mem
%   map  - the value of the map annotation
%
% GETMAPANNOTATION accepts the following parameter-value pairs
%   'iegm'    {:} | integer | array
%        - The electrogram point(s) for which the map annotation is
%        required
%
% GETMAPANNOTATION Returns the value(s) of the map annotation.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% woi = getReferenceAnnotation(userdata)
% woi = getReferenceAnnotation(userdata, 'iegm', 1)
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

ref = [];
if isstr(iEgm) && strcmpi(iEgm, ':')
    map = userdata.electric.annotations.mapAnnot;
else
    for i = 1:numel(iEgm)
        map(i,1) = userdata.electric.annotations.mapAnnot(iEgm(i),1);
    end
end

end