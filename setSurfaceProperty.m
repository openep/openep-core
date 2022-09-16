function newUserdata = setSurfaceProperty(userdata, varargin)
% SETSURFACEPROPERTY is used to set surface property values within OpenEP
% datasets
%
% Usage:
%   newUserdata = setSurfaceProperty(userdata)
% Where:
%   userdata - an OpenEP data structure
%   newUserdata - a new OpenEP data structure
%
% SETSURFACEPROPERTY accepts the following parameter-value pairs
%   'name'     {[]}|string
%   'map'      {[]}|array
%   'propSettings'  String containing the user-defined settings
%   'definedOn' {[]} | 'points' | 'elements'
%
% SETSURFACEPROPERTY Detailed description goes here
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

% get the parameter value pairs
nStandardArgs = 1;
name = [];
map = [];
propSettings = [];
definedOn = [];
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'name'
                name = varargin{i+1};
            case 'map'
                map = varargin{i+1};
            case 'propSettings'
                propSettings = varargin{i+1};
            case 'definedOn'
                definedOn = varargin{i+1};
        end
    end
end

% check if a property with that name exists and give a warning if so
if isfield(userdata.surface, 'surfaceProperties')
    allProperties = userdata.surface.surfaceProperties;
    if ~isempty(allProperties)
        for i = 1:numel(allProperties)
            if strcmpi(allProperties{i}.name, name)
                warning(['OPENEP/SURSURFACEPROPERTY: A surface property with name, ' name ' already exists. Data will be overwritten.']);
            end
        end
    end
end

% construct the property
if isempty(name)
    newUserdata = userdata;
    warning('OPENEP/SETSURFACEPROPERTY: No property data was provided');
else
    prop.name = name;
    prop.map = map;
    prop.propSettings = propSettings;
    prop.definedOn = definedOn;

    % store the property
    newUserdata = userdata;
    if ~isfield(userdata.surface, 'surfaceProperties')
        newUserdata.surface.surfaceProperties = {};
    end
    newUserdata.surface.surfaceProperties{end+1} = prop;
end

end