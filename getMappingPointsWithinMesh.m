function iPoint = getMappingPointsWithinMesh( userdata, varargin )
% GETMAPPINGPOINTSWITHINMESH Returns the indices of the mapping points
% which are located internal to the mesh
%
% Usage:
%   iPoint = getMappingPointsWithinMesh( userdata )
% Where:
%   userdata  - see importcarto_mem
%   iPoint  - logical array list of valid points; indexes into userdata.electric
%
% GETMAPPINGPOINTSWITHINMESH accepts the following parameter-value pairs
%   'tol' 0.1 | double
%       - The distance threshold within which points are considered to be
%         internal or external to the triangulation
%   'plot'     {false}|true
%       - Specify whether to plot the results
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

[iPoint, ~] = pointStatus( userdata, varargin{:} );

end