function op = getEgmTriangulation(userdata, varargin)
% GETEGMTRIANGULATION returns a triangulation of electrogram recording locations
%
% Usage:
%   tri = getEgmTriangulation( userdata )
% Where:
%   userdata - an OpenEP data structure
%   op - the triangulation, as a triangulation object, TriRep object or
%        structure, depending on the parameter 'type' value or as default a
%        triangulation
%
% GETEGMTRIANGULATION accepts the following parameter-value pairs
%   'type'     {'triangulation'} | 'trirep' | 'struct'
%
% GETEGMTRIANGULATION returns a triangulation representing the boundary of
% all the maping points. 
%
% Author: Steven Williams (2022)
% Modifications -
%
% See also BOUNDARY
%
% Info on Code Testing:
% ---------------------------------------------------------------
% load openep_dataset_1.mat
% tri = getEgmTriangulation(userdata); % returns a triangualtion object by default
% trisurf(tri, 'facecolor', [.4 .4 .4], 'edgecolor', 'k');
% axis equal vis3d; 
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1;
type = 'triangulation';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'type'
                type = varargin{i+1};
        end
    end
end

shrinkFactor = 1;
X = getElectrogramX(userdata);
tri = boundary(getElectrogramX(userdata), shrinkFactor);

switch type
    case 'triangulation'
        op = triangulation(tri, X);

    case 'trirep'
        op = TriRep(tri, X);

    case 'struct'
        op.X = X;
        op.Triangulation = tri;

end