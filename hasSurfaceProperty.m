function tf = hasSurfaceProperty(userdata, name)
% HASSURFACEPROPERTY returns whether the specified property exists
%
% Usage:
%   tf = hasSurfaceProperty(userdata, name)
% Where:
%   userdata - an OpenEP data structure
%   name     - the name of the required property
%
% HASSURFACEPROPERTY Detailed description goes here
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

tf = false;
if ~isfield(userdata.surface, 'surfaceProperties')
    warning('OPENEP/GETSURFACEPROPERTY: There are no surface properties defined in userdata');
else
    allProperties = userdata.surface.surfaceProperties;
    if ~isempty(allProperties)
        for i = 1:numel(allProperties)
            if strcmpi(allProperties{i}.name, name)
                tf = true;
            end
        end
    end
end

end