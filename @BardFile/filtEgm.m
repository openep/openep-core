function [fe, e] = filtEgm(hB, iTime, iChannel, varargin)
% @BARDFILE/FILTEGM Returns the filtered electrogram(s).
% Usage:
%   [fe, e] = filtEgm(hB, iTime, iChannel)
%   [fe, e] = filtEgm(hB, iTime, iChannel, 'grad')
% Where:
%   fe is the filtered electrogram
%   e is the raw electrogram (zero boundary conditions at limits of iTime)
%   'grad' - makes fe and e correspond to the gradient of the electrogram
% Author: Nick Linton (2012)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

isGrad = false;
if nargin==4
    if strstartcmpi('grad',varargin{1})
        isGrad = true;
    end
end

e = egm(hB, iTime, iChannel);
if isGrad
    e = e(3:end,:) - e(1:(end-2),:);
    e = e*(hB.SampleRate/2);
end
fe = zeros(size(e));

for i = 1:size(fe,2)
    fe(:,i) = filtfilt(hB.DFilt.sosMatrix,hB.PrivateDFilt.ScaleValues,e(:,i));
end

if isGrad
    pad = NaN(1,size(e,2));
    fe = [pad ; fe ; pad];
    e = [pad ; e ; pad];
    %warning('not sure this works')
end
