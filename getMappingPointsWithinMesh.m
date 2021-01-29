function iPoint = getMappingPointsWithinMesh( userdata, varargin )
% GETMAPPINGPOINTSWITHINMESH Returns the indices of the mapping points
% which are located internal to the mesh
%
% Usage:
%   iPoint = getMappingPointsWithinMesh( userdata )
% Where:
%   userdata  - see importcarto_mem
%   iPoint  - logical array list of valid points; indexes into
%             `userdata.electric`.
%
% GETMAPPINGPOINTSWITHINMESH accepts the following parameter-value pairs
%   'tol' 0.1 | double
%       - The distance threshold within which points are considered to be
%         internal or external to the triangulation
%   'plot'     {false}|true
%       - Specify whether to plot the results
%   'mask' {[]} | {Triangulation}
%       - A Triangulation object describing a supplementary three
%         dimensional mask to be used for filtering points
%
% GETMAPPINGPOINTSWITHINMESH Returns the indices of the mapping points
% which are located internal to the mesh. The OpenEP function pointStatus.m
% is used to identify these points. The parameter value pairs are passed
% directly onto pointStatus.m
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% See also GETMAPPINGPOINTSWITHINWOI, POINTSTATUS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% iPoint = getMappingPointsWithinMesh( userdata );
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
plot = false;
tol = 0.1;
mask = [];

if nargin > nStandardArgs && ~isempty(varargin{1})
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'plot'
                plot = varargin{i+1};
            case 'tol'
                tol = varargin{i+1};
            case 'mask'
                mask = varargin{i+1};
        end
    end
end

[iPoint, ~] = pointStatus( userdata, 'tol', tol, 'plot', plot);

if ~isempty(mask)
    % we need to find the points which are internal to the mask
    FV.vertices = mask.Points;
    FV.faces = mask.ConnectivityList;
    QPTS = userdata.electric.egmX(iPoint,:);
    iPoint = inpolyhedron(FV, QPTS, 'tol', tol);
end

end