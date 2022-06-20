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

if isa(userdata.surface.triRep, 'TriRep')
    vertices = double(userdata.surface.triRep.X);
elseif isa(userdata.surface.triRep, 'triangulation')
    vertices = double(userdata.surface.triRep.Points);
elseif isa(userdata.surface.triRep, 'struct')
    vertices = double(userdata.surface.triRep.X);
else
    error('OPENEP/getVertices: userdata.surface.TriRep should be a TriRep, structure or a triangulation object')
end
[~, isVertUsed] = repack(getMesh(userdata));

if used
    vertices = vertices(isVertUsed,:);
    isVertUsed = true(size(vertices));
end