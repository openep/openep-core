function activationTime = getActivationTime( userdata, varargin )
% GETACTIVATIONTIME Returns the activation times
%
% Usage:
%   activationTime = getActivationTime( userdata )
% Where:
%   userdata  - see importcarto_mem
%   activationTime  - the activation time for each electrogram
%
% GETACTIVATIONTIME accepts the following parameter-value pairs
%   'method'     {'annotation'} | 'nleo'
%   'units'      {'time'}       | 'samples'
%
% GETACTIVATIONTIME Returns activation times for each bipolar
% electrogram in userdata. The method depends on `method` as follows.
% (1) (... 'method', 'annotation')
%   This method returns the local activation times, relative to the
%   reference annotation, applied by the clinical system
% (2) (... 'method', 'nleo')
%   First, we iterate through all the bipolar mapping points in
%   userdata.electric. For each mapping point we get the electrogram within
%   the window of interest. We apply the non-linear energy operator to this
%   electrogram. We identify the earliest activation of the energy operator
%   and assign that to activationTime.
% Both methods use the sampling frequency stored in `userdata.electric.sampleFrequency
% and assume that sampling frequencey is 1000Hz if not specified. If the
% input parameter `'units'` is set to `'time'` (the default) then times are
% returned in seconds. If `'units'` is set to 'samples' then samples are
% returned.
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% activationTime = getActivationTime(userdata)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1;
method = 'annotation';
units = 'time';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'method'
                method = lower(varargin{i+1});
            case 'units'
                units = lower(varargin{i+1});
        end
    end
end

% Get the sampling frequency
sF = getSampleFrequency(userdata, 'default', 1000);

% Calculate the activation times
switch method
    % Using annotation data
    case 'annotation'
        activationTime = ( ...
            getMapAnnotation(userdata) -  ...
            getReferenceAnnotation(userdata) ...
            );
        
        % Using non-linear energy operator
    case 'nleo'
        iPoint = getMappingPointsWithinWoI(userdata);
        activationTime = [];
        for i = 1:numel(userdata.electric.names)
            if iPoint(i)
                woi = getWindowOfInterest(userdata, 'iegm', i);
                reference = getReferenceAnnotation(userdata, 'iegm', i);
                
                electrogram = getEgmsAtPoints(userdata, 'iegm', i, 'egmtype', 'bip', 'reference', 'off');
                e = electrogram{1}(reference+woi(1):reference+woi(2));
                
                [~,~,act] = nleo(e', 'adaptive', false, 'threshold', 5E-3, 'sampleFreq', sF);
                
                A = find(act);
                if ~isempty(A)
                    activationTime(i,1) = A(1); %#ok<*AGROW>
                else
                    warning(['OPENEP/EGMDURATION: No activation identified for point ' num2str(i)]);
                    activationTime(i,1) = NaN;
                end
            end
        end
end

% convert to time units if needed
if strcmpi(units, 'time')
    activationTime = activationTime / sF;
end
end