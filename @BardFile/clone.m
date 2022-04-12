function newObj = clone(oldObj)
% CLONE makes a new and separate copy of an object.
% Usage:
%   newObj = clone(oldObj)
%
% Author: Nick Linton (2011)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------
                        
    newObj = eval(class(oldObj));

    meta = ?BardFile;

    ndp = findobj([meta.Properties{:}],'Dependent',false);
    for idx = 1:length(ndp)
        newObj.(ndp(idx).Name) = oldObj.(ndp(idx).Name);
    end
end