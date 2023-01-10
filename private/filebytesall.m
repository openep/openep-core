function bytes = filebytesall(fID)
% FBYTES2END Gets number of bytes from start to end of file.
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

    oldPos = ftell(fID);
    if oldPos == (-1)
        error('Position in file could not be determined.')
    end
    fseek(fID, 0,'eof');
    bytes = ftell(fID);
    
    fseek(fID, oldPos, 'bof');

    