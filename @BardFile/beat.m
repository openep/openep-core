function [e, precedingCL, startIndex] = beat(hB, nBeat, iChannel, varargin)
% @BARDFILE/BEAT Accesses a particular beat from the run.
% Usage:
%   [e, precedingCL] = beat(hP, nBeat, iChannel)
%   [fe, precedingCL] = beat(hP, nBeat, iChannel, 'filter')
%   [{e,fe}, precedingCL] = beat(hP, nBeat, iChannel, 'both')
% Author: Nick Linton (2012)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------


stim = hB.StimIndices;

if isempty(stim)
    e = [];
    precedingCL = [];
    return
end

stim = [stim(:) ; stim(end)+round(hB.SampleRate*0.25)];

if numel(nBeat)>1
    error('@BARDFILE/BEAT: not coded to return more than one beat at a time.')
end

if nBeat==1
    precedingCL = NaN;
else
    precedingCL = (stim(nBeat) - stim(nBeat-1)) / hB.SampleRate;
end

if nargin == 4 && strcmpi(varargin{1}, 'filter')
    e = filtEgm(hB, stim(nBeat):stim(nBeat+1), iChannel);
elseif nargin == 4 && strcmpi(varargin{1}, 'both')
    [eFilt, eOrig] = filtEgm(hB, stim(nBeat):stim(nBeat+1), iChannel );
    e = {eOrig, eFilt};
else
    e = egm(hB, stim(nBeat):stim(nBeat+1), iChannel);
end

startIndex = stim(nBeat);