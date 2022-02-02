function getUACBoundaryConditions( userdata, varargin )
% GETUACBOUNDARYCONDITIONS Returns the boundary conditions for UAC
% calculation
%
% Usage:
%   getUACBoundaryConditions( userdata, varargin )
% Where:
%   userdata   - an OpenEP data structure
%
% GETUACBOUNDARYCONDITIONS accepts the following parameter-value pairs
%   'plot'     {false}|true
%
% GETUACBOUNDARYCONDITIONS identifies the boundary conditions for UAC calculation.
%
% Author: Steven Williams (2022) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% See also FREEBOUNDARYPOINTS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

figure
hold on
hSurf = drawMap(userdata, 'type', 'none');

[FF, l, a, tr] = getAnatomicalStructures( userdata, 'plot', false);

drawFreeBoundaries(FF, getMesh(userdata));

x = input('Which boundary is the SVC? ');
x = input('Which boundary is the IVC? ');
x = input('Which boundary is the CS? ');
x = input('Which boundary is the TV? ');

surfaceMesh = getMesh(userdata, 'type', 'triangulation', 'repack', true);
pp = PointPicker(surfaceMesh, hSurf, gca());
disp('select a point near the superior vena cava and a point near the inferior vena cava');
pause();

pointIndices = pp.PointIndices;

% Initialise the geodesic library and algorithm
initialiseGeodesic()
mesh = geodesic_new_mesh(surfaceMesh.Points, surfaceMesh.ConnectivityList);
algorithm = geodesic_new_algorithm(mesh, 'exact');

source_points{1} = geodesic_create_surface_point('vertex', pointIndices(1), surfaceMesh.Points(pointIndices(1),:));
geodesic_propagate(algorithm, source_points);
destination = geodesic_create_surface_point('vertex', pointIndices(2), surfaceMesh.Points(pointIndices(2),:));
path = geodesic_trace_back(algorithm, destination);


[x,y,z] = extract_coordinates_from_path(path);




plot3( x*1.001 ...
    ,y*1.001 ...
    ,z*1.001 ...
    ,'k-s' ...
    , 'LineWidth', 2 ...
    , 'markersize', 10 ...
    );


geodesic_delete();

end