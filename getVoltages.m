function [voltages] = getVoltages( userdata, varargin )
% GETVOLTAGES Returns the voltages associated with each mapping point
%
% Usage:
%   [voltages] = getVoltages( userdata )
% Where:
%   userdata  - see importcarto_mem
%   voltages - all the voltages
%
% GETVOLTAGES accepts the following parameter-value pairs
%   'type'     {'bip'} | 'uni'
%
% GETVOLTAGES Returns the voltage values associated with each mapping
% point. Bipolar or unipolar voltages can be specified by setting the
% `type` property.
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% voltages = getVoltages( userdata, 'type', 'bip' )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
type  = 'bip';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'type'
                type = varargin{i+1};
        end
    end
end

switch lower(type)
    case 'bip'
        voltages = userdata.electric.voltages.bipolar;
    case 'uni'
        voltages = userdata.electric.voltages.unipolar;
end

end