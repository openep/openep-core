function V = calculatePeak2PeakVoltage( egms, refAnnot, woi )
% CALCULATEPEAK2PEAKVOLTAGE Calcualtes peak to peak voltages
%
% Usage:
%   V = calculatePeak2PeakVoltage( egms, refAnnot, woi )
% Where:
%   egms      - the electrograms
%   refAnnot  - the reference annotation, in samples
%   woi       - the window of interest, relative to the refAnnot, in samples
%   V         - the output voltages
%
% CALCULATEPEAK2PEAKVOLTAGE does not accept parameter-value pairs
%
% CALCULATEPEAK2PEAKVOLTAGE Detailed description goes here
%
% Author: Steven Williams (2025)
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

nEgm = size(egms,1);
V = NaN(nEgm,1);
for iEgm = 1:nEgm

    sampleRange = refAnnot(iEgm)+woi(iEgm,1):refAnnot(iEgm)+woi(iEgm,2);
    thisEgm = egms(iEgm, sampleRange);
    V(iEgm) = max(thisEgm) - min(thisEgm);

end
end
