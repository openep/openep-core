function [vertices, isVertUsed] = getVertices( userdata, varargin )
% GETVERTICES Returns the vertices referenced by userdata
%
% Usage:
%   [vertices, vertsref] = getVertices( userdata )
% Where:
%   userdata  - see importcarto_mem
%   vertices - all the vertices
%   isVertUsed - whether the vertex is referenced by the triangulation
%
% GETVERTICES accepts the following parameter-value pairs
%   'used'     {false} | true
%
% GETVERTICES Returns the vertices referenced by userdata. If the property
% `used` is set to `true` then only vertices in the triangulation that are
% used by the triangulation are returned.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%   SW 09-04-2021: add routine to only return used vertices
%
% Info on Code Testing:
% ---------------------------------------------------------------
% [vertices, isVertUsed] = getVertices( userdata )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
used  = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'used'
                used = varargin{i+1};
        end
    end
end

% Check the type of surface data stored in the OpenEP datastructure and
% accordingly access the points
if isa(userdata.surface.triRep, 'TriRep')
    FV.vert = userdata.surface.triRep.X;
    FV.faces = userdata.surface.triRep.Triangulation;
elseif isa(userdata.surface.triRep, 'triangulation')
    FV.vert = userdata.surface.triRep.Points;
    FV.faces = userdata.surface.triRep.ConnectivityList;
elseif isa(userdata.surface.triRep, 'struct')
    if ~isfield(userdata.surface.triRep, 'X')
        error('OPENEP:invalidData', 'OPENEP/GETMESH: invalid data. userdata.surface.triRep must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
    if ~isfield(userdata.surface.triRep, 'Triangulation')
        error('OPENEP:invalidData', 'OPENEP/GETMESH: invalid data. userdata.surface.triRep must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
    FV.vert = userdata.surface.triRep.X;
    FV.faces = userdata.surface.triRep.Triangulation;
else
    error('OPENEP:invalidData', 'OPENEP/GETMESH: userdata.surface.triRep must be one of: TriRep, triangulation or struct');
end

vertices = FV.vert;
[~, isVertUsed] = repack(getMesh(userdata));

if used
    vertices = vertices(isVertUsed,:);
    isVertUsed = true(size(vertices));
end

end