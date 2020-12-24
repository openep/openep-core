function trf = strstartcmpi(startpattern, str )
% STRSTARTCMPI Compares the begining STR with the other STARTPATTERN.
% Usage:
%   trf = strstartcmp(startpattern, str)
% Where:
%   str and startpattern are strings, str can be a cellstr
%
% STRSTARTCMPI Is similar to STRSTARTCMP excspt that the result is
% independent of case.

% Author: Nick Linton (2009)
% Modifications - 
% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

trf = strstartcmp(lower(startpattern), lower(str));



