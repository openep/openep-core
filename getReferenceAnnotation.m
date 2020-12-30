function ref = getReferenceAnnotation(userdata, varargin)
% GETREFERENCEANNOTATION Returns the value of the reference annotation 
%
% Usage:
%   ref = getReferenceAnnotation( userdata )
% Where:
%   userdata  - see importcarto_mem
%   ref  - the value of the reference annotation
%
% GETREFERENCEANNOTATION accepts the following parameter-value pairs
%   'iegm'    {:} | integer | array
%        - The electrogram point(s) for which the reference annotation is
%        required
%
% GETREFERENCEANNOTATION Returns the value fo the reference annotation.
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
    ref = userdata.electric.annotations.referenceAnnot;
else
    for i = 1:numel(iEgm)
        ref(i,1) = userdata.electric.annotations.referenceAnnot(iEgm(i),1);
    end
end

end