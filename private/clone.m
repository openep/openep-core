function newObj = clone(oldObj)
% NOTE: TO USE THIS, YOU MAY NEED TO COPY IT TO THE OBJECT'S HOME FOLDER
%
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
    
    warning('this code does not copy all properties from super/sub classes')
    newObj = eval(class(oldObj));

    meta = eval(['?',class(oldObj)]);

    ndp = findobj([meta.Properties{:}],'Dependent',false);
    for idx = 1:length(ndp)
        newObj.(ndp(idx).Name) = oldObj.(ndp(idx).Name);
    end
end