function prop = getSurfaceProperty(userdata, name)
% GETSURFACEPROPERTY returns the required surface property
%
% Usage:
%   prop = getSurfaceProperty(userdata, name)
% Where:
%   userdata - an OpenEP data structure
%   name     - the name of the required property
%   prop     - the surface property
%
% GETSURFACEPROPERTY Detailed description goes here
%
% Author: Steven Williams (2022)
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

% check if a property with that name exists and give a warning if so
if ~isfield(userdata.surface, 'surfaceProperties')
    warning('OPENEP/GETSURFACEPROPERTY: There are no surface properties defined in userdata');
else
    allProperties = userdata.surface.surfaceProperties;
    if ~isempty(allProperties)
        for i = 1:numel(allProperties)
            if strcmpi(allProperties{i}.name, name)
                prop = allProperties{i};
            end
        end
    end
end

end