function [nl, nlfilt, act] = nleo(x, varargin)
% NLEO returns the non-linear energy operator of a signal.
% Usage:
%   nl = nleo(x)
%   [nl, nlfilt, act] = nleo(x)
%
% Where:
%   x - is the input signal
%   nleo - is the output (nleo) signal
%   nleofilt - is the output filtered signal
%   act - is binary activation, defined by an adaptive algorithm
%
% NLEO uses the non-linear energy operator defined by Kaiser:
%   E_j = x_j^2 - (x_j+1)(x_j-1)
%
% Author: Steven Williams (2012)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% testcode
% ---------------------------------------------------------------
%
% Use warnings with the following format
%    disp('NLEO: warning.')
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

adaptive = true;
sampleFreq = 2000;

% calculate the parameters
xsq = x.^2;
x_jp = [0; x(1:end-1)];
x_jm = [x(2:end); 0];

% top and tale the signals
xsq = xsq(2:end-1);
x_jp = x_jp(2:end-1);
x_jm = x_jm(2:end-1);

% calculate the nleo
nl = xsq - (x_jp .* x_jm);

% add to the top and tail
nl = [0; nl; 0];

% low pass filter the nleo
cutOff = 24;                 %cutoff frequency
Wn = cutOff/(sampleFreq/2);  %normalised cutoff frequency
n = sampleFreq/200;          %filter order
b = fir1(n, Wn, 'low');      %design a lowpass filter
nlfilt = filtfilt(b, 1, nl); %zero-phase digital filtering

N = 100;
Fp = 24; % 24Hz passband frequency
Fs = sampleFreq; % 2000Hz sample frequency


% STATIC determine the activation
if nargin == 2
    t = varargin{2};
else
    t = 1E-4;
end

act = nlfilt;
act(act>t) = 1;
act(act<=t) = 0;

% ADAPTIVE determine the activation
if adaptive
    
    % set the window length
    window = 200; %ms
    windowLength = window/1000 * sampleFreq;
    stepLength = 10 * sampleFreq/1000;
    
    % initialise variables
    i = 1;
    iStart = ((i-1) * stepLength) + 1;
    iEnd = iStart + windowLength;
    
    t = nan(windowLength/stepLength,length(x));
    
    %weighting factor
    k = 2;
    
    % loop through the windows and assign the threshold
    
    while iEnd<length(x)
        t(i,iStart:iEnd) = k * std(nlfilt(iStart:iEnd));
        i = i+1;
        iStart = iStart + stepLength;
        iEnd = iEnd + windowLength;
    end
        
    % take the smallest of each thresh value
    threshold = min(t)';
    
    % top and tail the signal
    %threshold = threshold(2:end-1);
    
    act = nlfilt;
    act(act>=threshold) = 1;
    act(act<=threshold) = 0;
    
end
end
