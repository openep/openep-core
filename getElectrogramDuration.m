function egmDuration = getElectrogramDuration( userdata )
% GETELECTROGRAMDURATION Returns the electrogram durations
%
% Usage:
%   egmDuration = getElectrogramDuration( userdata )
% Where:
%   userdata  - see importcarto_mem
%   egmDuration  - the duration of activation for each electrogram
%
% GETELECTROGRAMDURATION Calculates electrogram durations for each bipolar
% electrogram in userdata. The method is as follows. First, we iterate
% through all the bipolar mapping points in userdata.electric. For each
% mapping point we get the electrogram within the window of interest. We
% apply the non-linear energy operator to this electrogram. We identify the
% earliest and latest activation of the energy operator. We calculate the
% difference and assign that to egmDuration. We also remove any values that
% are not within +/- 2 standard deviations of the mean. This is an
% arbritary cut off but some filtering is necessary to avoid identifying
% noise as continuous activation, for example.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% egmDuration = getElectrogramDuration(userdata)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

iPoint = getMappingPointsWithinWoI(userdata);

egmDuration = [];
for i = 1:numel(userdata.electric.names)
    if iPoint(i)
        woi = getWindowOfInterest(userdata, 'iegm', i);
        reference = getReferenceAnnotation(userdata, 'iegm', i);
        
        electrogram = getEgmsAtPoints(userdata, 'iegm', i, 'egmtype', 'bip', 'reference', 'off');
        e = electrogram{1}(reference+woi(1):reference+woi(2));
        
        [~,~,act] = nleo(e', 'adaptive', false);
        
        A = find(act);
        if ~isempty(A)
            egmDuration(i,1) = A(end) - A(1);
        else
            warning(['OPENEP/EGMDURATION: No activation identified for point ' num2str(i)]);
            egmDuration(i,1) = NaN;
        end  
    end
end

% remove any values that do not fall wthin ±2 standard deviations of the
% mean
m = nanmean(egmDuration);
s = nanstd(egmDuration);
toRemove = false(size(egmDuration));
toRemove(egmDuration<m-2*s) = 1;
toRemove(egmDuration>m+2*s) = 1;
egmDuration(toRemove) = NaN;

end