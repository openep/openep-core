function meanVoltage = getMeanVoltage( userdata, varargin )
% GETMEANVOLTAGE Returns the mean voltage
%
% Usage:
%   meanVoltage = getMeanVoltage( userdata, varargin )
% Where:
%   userdata  - see importcarto_mem
%   meanVoltage  - the mean chamber voltage (mV)
%
% GETMEANVOLTAGE accepts the following parameter-value pairs
%   'method'    {'map'} | 'egm'
%   'type'      {'bip'} | 'uni'
%
% GETMEANVOLTAGE Returns the mean voltage of a chamber. By default, the
% mean bipolar voltage is calculated using the interpolated mapping data
% from the clinical mapping system (stored in userdata.surface.act_bip). If
% 'method' is set to 'egm' then the bipolar voltage is first interpolated 
% from the bipolar electrogram data (stored in userdata.electric). If 
% 'type' is set to 'uni' then unipolar voltages are returned.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
method = 'map';
type = 'bip';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'type'
                type = varargin{i+1};
        end
    end
end

switch method
    case 'map'
        meanVoltage = nanmean(userdata.surface.act_bip(:,2));
    case 'egm'
        interpData = generateInterpData(userdata, 'bip-map');
        meanVoltage = nanmean(interpData);
end

