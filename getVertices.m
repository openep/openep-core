function [vertices, isVertUsed] = getVertices( userdata )
% GETVERTICES Returns the vertices referenced by userdata
%
% Usage:
%   [vertices, vertsref] = getVertices( userdata )
% Where:
%   userdata  - see importcarto_mem
%   vertices - all the vertices
%   isVertUsed - whether the vertex is referenced by the triangulation
%
% GETVERTICES Returns the vertices referenced by userdata
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% [vertices, isVertUsed] = getVertices( userdata )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

tr = getMesh(userdata, 'type', 'triangulation');
vertices = tr.X;
[~, isVertUsed] = repack(tr);

end
