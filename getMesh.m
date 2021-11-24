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
repack = 'false';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'type'
                type = varargin{i+1};
            case 'repack'
                repack = varargin{i+1};
        end
    end
end
if ~any(strcmpi({'trirep' 'triangulation'}, type))
    error(['OPENEP/GETMESH: Value: ' type ' for parameter: type not recognised']);
end

FV.vert = userdata.surface.triRep.X;
FV.faces = userdata.surface.triRep.Triangulation;

switch lower(type)
    case 'trirep'
        if isa(userdata.surface.triRep, 'TriRep')
            tr = userdata.surface.triRep;
            return;
        else
            tr = TriRep(FV.faces, FV.vert(:,1), FV.vert(:,2), FV.vert(:,3));
            return;
        end
    case 'triangulation'
        if isa(userdata.surface.triRep, 'triangulation')
            tr = userdata.surface.triRep;
            return;
        else
            tr = triangulation(FV.faces, FV.vert(:,1), FV.vert(:,2), FV.vert(:,3));
            return;
        end
end

if repack
    tr = repack(tr);
end
