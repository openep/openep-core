function newUserdata = setMesh(userdata, tNew)
% SETMESH Stores the new mesh tNew in userdata
%
% Usage:
%   newUserdata = setMesh(userdata, tNew)
% Where:
%   userdata - an OpenEP dataset, see https://openep.io/data/
%   tNew - the new mesh. One of;
%           a Matlab TriRep object
%           a Matlab triangulation object
%           a structure with the fields .X and .Triangulation
%
% SETMESH Returns a copy of userdata with the mesh changed for the new
% mesh. Note that the other fields in userdata.surface representing the
% datatypes will therfore become invalid unless further adjustments are
% made to these fields outwith this function.
%
% See also REPACKUSERDATA.
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
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

% Check the input data
if isa(tNew, 'TriRep')
    type = 'trirep';
elseif isa(tNew, 'triangulation')
    type = 'triangulation';
elseif isa(tNew, 'struct')
    type = 'struct';
    if ~isfield(userdata.surface.triRep, 'X')
        error('OPENEP/SETMESH: invalid input data. tNew must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
    if ~isfield(userdata.surface.triRep, 'Triangulation')
        error('OPENEP/SETMESH: invalid intput data. tNew must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
end



FV.vert = userdata.surface.triRep.X;
FV.faces = userdata.surface.triRep.Triangulation;

switch lower(type)
    case 'trirep'
        if isa(userdata.surface.triRep, 'TriRep')
            tr = userdata.surface.triRep;
        else
            tr = TriRep(FV.faces, FV.vert(:,1), FV.vert(:,2), FV.vert(:,3));
        end
    case 'triangulation'
        if isa(userdata.surface.triRep, 'triangulation')
            tr = userdata.surface.triRep;
        else
            tr = triangulation(FV.faces, FV.vert(:,1), FV.vert(:,2), FV.vert(:,3));
        end
end

if dorepack
    tr = repack(tr);
end
