function cvHistogram( userdata, varargin )
% CVHISTOGRAM Draws a conduction velocity histogram
%
% Usage:
%   cvHistogram( userdata )
% Where:
%   userdata  - see importcarto_mem
%
% CVHISTOGRAM accepts the following parameter-value pairs
%   'limits'    {[0 5]} | array
%   'binwidth'  {0.1} | double
%
% CVHISTOGRAM displays a histogram of conduction velocities. Limits are set
% to exclude non-physiological conduction velocities
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% cvHistogram( userdata )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
limits = [0 5];
binwidth = 0.1;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'limits'
                limits = varargin{i+1};
            case 'binwidth'
                binwidth = varargin{i+1};
        end
    end
end
% TODO: check format of input values for each parameter are correct

cvdata = getConductionVelocity(userdata);
cvdata(cvdata<limits(1)) = [];
cvdata(cvdata>limits(2)) = [];

histogram(cvdata ...
    , 'normalization', 'probability' ...
    , 'binwidth', binwidth ...
    , 'facecolor', 'k' ...
    , 'facealpha', 1 ...
    );
set(gcf, 'color', 'w');
xlabel('Conduction Velocity (m/s)');
ylabel('Frequency')
set(gca, 'fontsize', 16)
end
