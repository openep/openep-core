function sampleFrequency = getSampleFrequency( userdata, varargin )
% GETSAMPLEFREQUENCY Returns the electrogram sampling frequency
%
% Usage:
%   sampleFrequency = getSampleFrequecny( userdata )
% Where:
%   userdata  - see importcarto_mem
%   sampleFrequency  - the electrogram sample frequency
%
% GETSAMPLEFREQUENCY accepts the following parameter-value pairs
%   'default'     {NaN} | 1000
%
% GETSAMPLEFREQUENCY Returns the sample frequency of the electrograms. If
% no sample frequency is specified, then NaN is returned. If the parameter
% `default` has been set, then if sample frequency data is not found then 
% the value of `default` is returned as the sample frequency. 
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% sampleFrequency = getSampleFrequency(userdata)
% sampleFrequency = getSampleFrequency(userdata, 'default', 1000);
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------
  
nStandardArgs = 1; % UPDATE VALUE
default = NaN;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'default'
                default = varargin{i+1};
        end
    end
end

if isfield(userdata.electric, 'sampleFrequency')
    if isnumeric(userdata.electric.sampleFrequency)
        sampleFrequency = userdata.electric.sampleFrequency;
        return;
    else
        error('OPENEP/GETSAMPLEFREQUENCY: Unrecognised format of sampleFrequency in dataset');
    end
else
    if isnan(default)
        sampleFrequency = NaN;
    elseif isnumeric(default)
        sampleFrequency = default;
    else
        error('OPENEP/GETSAMPLEFREQUENCY: Unrecognised value for parameter `default`');
    end
end