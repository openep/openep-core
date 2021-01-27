function cmdout = du( argument )
% DU is a wrapper function for the Unix file system tool du
%
% Usage:
%   cmdout = du(argument)
% Where:
%   argument - is the input provided to do
%   cmdout - is the information returned by du in a string
%
% DU Detailed description goes here
%
% Author: Steven Williams (2020)
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

%create the system call
arg = ['du ' argument];

%system call
[status, cmdout] = system(arg);

%remove the carriage return from the end of cmdout
cmdout(end) = [];

end