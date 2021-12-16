function tr = getMesh(userdata, varargin)
% GETMESH Returns the triangulation-based mesh from userdata
%
% Usage:
%   tr = getMesh(userdata)
% Where:
%   tr - a TriRep, or Triangulation, object
%
% GETMESH accepts the following parameter-value pairs
%   'type'     {'trirep'}|'triangulation'|'struct'
%       - Specifies whether to return the mesh as a TriRep object or as a
%       Triangulation object or as a struct
%   'limitToTriangulation' {'false'}|true
%       - Specifies whether to repack the triangulation
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
limitToTriangulation = false;
dorepack = false;

if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'type'
                type = varargin{i+1};
            case 'limittotriangulation'
                limitToTriangulation = varargin{i+1};
            case 'repack'
                dorepack = varargin{i+1};
        end
    end
end
if ~any(strcmpi({'trirep' 'triangulation' 'struct'}, type))
    error(['OPENEP/GETMESH: Value: ' type ' for parameter: type not recognised']);
end

% Check the type of surface data stored in the OpenEP datastructure and
% accordingly access the triangulation and the points
if isa(userdata.surface.triRep, 'TriRep')
    FV.vert = userdata.surface.triRep.X;
    FV.faces = userdata.surface.triRep.Triangulation;
elseif isa(userdata.surface.triRep, 'triangulation')
    FV.vert = userdata.surface.triRep.Points;
    FV.faces = userdata.surface.triRep.ConnectivityList;
elseif isa(userdata.surface.triRep, 'struct')
    if ~isfield(userdata.surface.triRep, 'X')
        error('OPENEP/GETMESH: invalid data. userdata.surface.triRep must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
    if ~isfield(userdata.surface.triRep, 'Triangulation')
        error('OPENEP/SETMESH: invalid data. userdata.surface.triRep must be one of: TriRep, triangulation or struct with fields .X and .Triangulation');
    end
    FV.vert = userdata.surface.triRep.X;
    FV.faces = userdata.surface.triRep.Triangulation;
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
    case 'struct'
        if isa(userdata.surface.triRep, 'struct')
            tr = userdata.surface.triRep;
        else
            tr.X = FV.vert;
            tr.Points = FV.faces;
        end
end

if limitToTriangulation
    tr = repack(tr);
end
