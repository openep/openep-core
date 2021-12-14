function tr = getMesh(userdata, varargin)
% GETMESH Returns the triangulation-based mesh from userdata
%
% Usage:
%   tr = getMesh(userdata)
% Where:
%   tr - a TriRep, or Triangulation, object
%
% GETMESH accepts the following parameter-value pairs
%   'type'     {'trirep'}|'triangulation'
%       - Specifies whether to return the mesh as a TriRep object or as a
%       Triangulation object
<<<<<<< HEAD
%   'limitToTriangulation' {'false'}|true
%       - Specifies whether to repack the triangulation
=======
%   'repack'    {false}|true
>>>>>>> origin/develop
%
% GETMESH Returns a face/vertex representation of the anatomical model. 
% Supported data types include istances of the Matlab objects Trirep and 
% Triangulation.
%
% Author: Steven Williams (2020) (Copyright)
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

nStandardArgs = 1; % UPDATE VALUE
type = 'trirep';
<<<<<<< HEAD
limitToTriangulation = false;
=======
dorepack = false;
>>>>>>> origin/develop
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'type'
                type = varargin{i+1};
<<<<<<< HEAD
            case 'limittotriangulation'
                limitToTriangulation - varargin{i+1};
=======
            case 'repack'
                dorepack = varargin{i+1};
>>>>>>> origin/develop
        end
    end
end
if ~any(strcmpi({'trirep' 'triangulation'}, type))
    error(['OPENEP/GETMESH: Value: ' type ' for parameter: type not recognised']);
end

<<<<<<< HEAD
=======
% Check the type of surface data stored in the OpenEP datastructure and
% accordingly access the triangulation and the points
>>>>>>> origin/develop
if isa(userdata.surface.triRep, 'TriRep')
    FV.vert = userdata.surface.triRep.X;
    FV.faces = userdata.surface.triRep.Triangulation;
elseif isa(userdata.surface.triRep, 'triangulation')
    FV.vert = userdata.surface.triRep.Points;
    FV.faces = userdata.surface.triRep.ConnectivityList;
<<<<<<< HEAD
else
    error('OPENEP/getMesh: userdata.surface.TriRep should be a TriRep or a triangulation object')
=======
elseif isa(userdata.surface.triRep, 'struct')
    if ~isfield(userdata.surface.triRep, 'X')
        error('OPENEP/GETMESH: invalid data. userdata.surface.triRep must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
    if ~isfield(userdata.surface.triRep, 'Triangulation')
        error('OPENEP/SETMESH: invalid data. userdata.surface.triRep must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
    FV.vert = userdata.surface.triRep.X;
    FV.faces = userdata.surface.triRep.Triangulation;
>>>>>>> origin/develop
end

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

<<<<<<< HEAD
if limitToTriangulation
    tr = repack(tr);
end

=======
if dorepack
    tr = repack(tr);
>>>>>>> origin/develop
end
