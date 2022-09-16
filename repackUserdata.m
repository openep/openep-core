function newUserdata = repackUserdata(userdata)
% REPACKUSERDATA Removes unreferenced points from the triangulation and all
% point-based datatypes within OpenEP
%
% Usage:
%   newUserdata = repackUserdata( userdata )
% Where:
%   userdata - the input OpenEP dataset, see https://openep.io/data/
%   newUserdata  - the output OpenEP dataset
%
% REPACKUSERDATA removes all non-referenced points and components of point
% based datasets. These are stored within userdata.surface
%
% See also SETMESH
%
% Author: Steven Williams (2021)
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

% get the new mesh and identify the vertices not referenced by the triangulation
[tNew, isVertUsed] = repack(getMesh(userdata));

% set the new mesh
newUserdata = setMesh(userdata, tNew);

% remove irrelevant data
newUserdata.surface.act_bip(~isVertUsed,:) = [];
newUserdata.surface.uni_imp_frc(~isVertUsed,:) = [];
newUserdata.surface.isVertexAtRim(~isVertUsed,:) = [];

end