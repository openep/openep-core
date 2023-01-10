function bytes = filebytes2end(fID)
% FBYTES2END Gets number of bytes from current position to end of file.
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

    bytes = filebytesall(fID)- ftell(fID);
    