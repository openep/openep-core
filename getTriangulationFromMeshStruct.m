function tr = getTriangulationFromMeshStruct( meshStruct, varargin )
% GETTRIANGULATIONFROMMESHSTRUCT Returns a Triangulation object from a
% structure describing a mesh
%
% Usage:
%   tr = getTriangulationFromMeshStruct( meshStruct )
% Where:
%   meshStruct - a structure containing .Pts and .Tri fields. Note that the
%                final column in .Tri contains a region label
%
% GETTRIANGULATIONFROMMESHSTRUCT accepts the following parameter-value pairs
%   'region' {':'} | int | int array
%       - The required region from the mesh structure
%   'scale'  {'mm'} | 'um'
%       - The units of length for the meshStruct
%   'repack' {true} | false
%       - Whether to repack the triangulation to remove unused vertices
%
% GETTRIANGULATIONFROMMESHSTRUCT Creates a Triangulation object from a
% meshStructure which contains .Pts and .Tri fields. Note that
% GETTRIANGULATIONFROMMESHSTRUCT supports region labelling. The final
% column in .Tri therefore should contain region labels as integers. The
% parameter `region` can be specified such that on the specified sub
% region(s) are added to the Triangulation object. By default, repack() is
% called on the Triangulation such that redundant points are removed from
% the mesh. This also changes the vertex indexing. If this behaviour is not
% required then specify the parameter `repack` to be `false`. Units of the
% Triangulation.Points are mm, by default, but this will only be correct of
% the input parameter `scale` is set; or if the input units are mm.
% Possible input scales are mm or um.
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% tr = getTriangulationFromMeshStruct( meshStruct, 'region', 11 )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
region = ':';
dorepack = true;
scale = 'mm';
type = 'Triangulation';
if nargin > nStandardArgs && ~isempty(varargin{1})
    for i = 1:2:nargin-1
        switch lower(varargin{i})
            case 'region'
                region = varargin{i+1};
            case 'repack'
                dorepack = varargin{i+1};
            case 'scale'
                scale = varargin{i+1};
            case 'type'
                type = varargin{i+1};
                
        end
    end
end

if strcmpi(region, ':')
    % then we want the whole mesh structure
    T = meshStruct.Tri(:,1:3);
else
    T = meshStruct.Tri(meshStruct.Tri(:,4)==region,1:3);
end

switch lower(scale)
    case 'mm'
        meshStruct.Pts = meshStruct.Pts;
    case 'um'
        meshStruct.Pts = meshStruct.Pts / 1000;
end

switch lower(type)
    case 'triangulation'
        tr = triangulation(T, meshStruct.Pts(:,1), meshStruct.Pts(:,2), meshStruct.Pts(:,3));
    case 'trirep'
        tr = TriRep(T, meshStruct.Pts(:,1), meshStruct.Pts(:,2), meshStruct.Pts(:,3));
end

if dorepack
    tr = repack(tr);
end