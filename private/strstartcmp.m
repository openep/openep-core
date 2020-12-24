function trf = strstartcmp(startpattern, str)
% STRSTARTCMP Compares the begnining STR with the other STARTPATTERN.
% Usage:
%   trf = strstartcmp(startpattern, str)
% Where:
%   str and startpattern are strings, str can be a cellstr


% Author: Nick Linton (2009)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

    if iscell(str)
        trf = false(size(str));
        for i = 1:numel(str)
            trf(i) = mystrstartcmp(str{i});
        end
    else
        trf = mystrstartcmp(str);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = mystrstartcmp(str)
    if length(str) >= length(startpattern)
        tf = strncmp(startpattern, str, length(startpattern));
    else
        tf = false;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end





