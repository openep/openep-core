function [ivcIndices, svcIndices, csIndices, tvPointSet, postPointSet] = getUACBoundaryConditions( userdata, varargin )
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

% Draw the surface
figure
hold on
hSurf = drawMap(userdata, 'type', 'none');

% Get the free boundaries
[FF, l, a, tr] = getAnatomicalStructures( userdata, 'plot', false);

% Draw the free boundaries
drawFreeBoundaries(FF, getMesh(userdata));

% Create a legend with numbers for each boundary
legendText{1} = '';
legendText{2} = '';
for i = 1:numel(FF)
legendText{i+2} = ['Boundary ' num2str(i)];
end
legend(legendText)

% User input to label the boundaries
iSVC = input('Which boundary is the SVC? ');
iIVC = input('Which boundary is the IVC? ');
iCS = input('Which boundary is the CS? ');
iTV = input('Which boundary is the TV? ');

% update the legend
legendText{iSVC+2} = 'SVC';
legendText{iIVC+2} = 'IVC';
legendText{iCS+2} = 'CS';
legendText{iTV+2} = 'TV';
legend(legendText);

% get the point sets for each of the boundaries
ivcIndices = unique(FF{iIVC});
svcIndices = unique(FF{iSVC});
csIndices = unique(FF{iCS});
tvIndices = unique(FF{iTV});


% select a point near the SVC and a point near the IVC
originalMesh = getMesh(userdata, 'type', 'triangulation');
surfaceMesh = getMesh(userdata, 'type', 'triangulation', 'limittotriangulation', true);
pp = PointPicker(originalMesh, hSurf, gca());
disp('select a point near the superior vena cava and a point near the inferior vena cava');
pause();

% get these points - the first should be the SVC
pointIndices = pp.PointIndices;

% find the closest point in the SVC boundary for the SVC point

% So
% pointIndices(1) is the SVC point, indexing into originalMesh.Points
% pointIndices(2) is the IVC point, indexing into originalMesh.Points

% find the SVC point on the ring
svcPointOnRing = findclosestvertex(originalMesh.Points(svcIndices,:), originalMesh.Points(pointIndices(1),:));
svcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(svcIndices(svcPointOnRing),:));
 
ivcPointOnRing = findclosestvertex(originalMesh.Points(ivcIndices,:), originalMesh.Points(pointIndices(2),:));
ivcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(ivcIndices(ivcPointOnRing),:));

plotTag(userdata, 'coord', surfaceMesh.Points(svcPointOnSurface,:), 'color', 'g');
plotTag(userdata, 'coord', surfaceMesh.Points(ivcPointOnSurface,:), 'color', 'g');


% Initialise the geodesic library and algorithm
initialiseGeodesic()
mesh = geodesic_new_mesh(surfaceMesh.Points, surfaceMesh.ConnectivityList);
algorithm = geodesic_new_algorithm(mesh, 'exact');

source_points{1} = geodesic_create_surface_point('vertex', svcPointOnSurface, surfaceMesh.Points(svcPointOnSurface,:));
geodesic_propagate(algorithm, source_points);
destination = geodesic_create_surface_point('vertex', ivcPointOnSurface, surfaceMesh.Points(ivcPointOnSurface,:));
path = geodesic_trace_back(algorithm, destination);


[x,y,z] = extract_coordinates_from_path(path);

plot3( x*1.001 ...
    ,y*1.001 ...
    ,z*1.001 ...
    ,'k-' ...
    , 'LineWidth', 2 ...
    , 'markersize', 10 ...
    );

geodesic_delete();

end