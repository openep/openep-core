function bytes = filebytesbetweenindices(fID)
% FILEBYTESBETWEENINDICES Gets number of bytes from current position to end of file.
%
% Usage:
%   bytes = fbytes2end(fileID)
%
% Where:
%   fileID is the ID of an open file

% Author: Nick Linton (2015)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

fseek(fID, pos1, 'bof');
bytes(1) = ftell(fID);
fseek(fID, pos2, 'bof');
bytes(2) = ftell(fID);

bytes = bytes(2) - bytes(1);
    