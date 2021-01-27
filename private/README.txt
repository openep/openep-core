This includes some functions for global interpolation of scattered data using RBFs.

It can be used to approximate conduction velocities.

Given some points at which local activation times 'LATs' are defined, 
use "RBFConductionVelocity" to approximate conduction velocities anywhere. 

Example - copy and paste below the line:
----------------------------------------------------------------------

load('data/example.mat')

% size(rovingUsed)=[d,n] where d is dimensions and n is number of points
% size(LATs)=[1,n] (a scalar field)
% size(verts)=[N,d] where N does not have to equal n
% In this example 'verts' contains the vertices of a surface triangulation from Precision

[interpLATs, d_interpLATs, u, n, speed] = RBFConductionVelocity(LATs, rovingUsed, verts');

% Plot some results:

% Plot interpolated LATs

subplot(1,2,1),trisurf(connects,verts(:,1),verts(:,2),verts(:,3),interpLATs','EdgeColor','none');
axis equal;
colorbar;
title('Interpolated LATs');

% Overlay unit-vectors of the conduction velocity:

% first project the velocities in tangent direction
projU=projectVertexVectors(verts,u',normals);

% then compute unit-vectors:
uv=createUnitVectors(projU);

hold on
quiver3(verts(:,1),verts(:,2),verts(:,3),uv(:,1),uv(:,2),uv(:,3),'color','k');
hold off


% Also plot approximated conduction speed:

subplot(1,2,2),trisurf(connects,verts(:,1),verts(:,2),verts(:,3),speed','EdgeColor','none');
axis equal;
colorbar; % units here are distance/time, in units defined by 'LATs' and 'verts'
title('Interpolated speed (distance/time)');

% Overlay the nearest-neighbour points on the mesh
hold on
plot3(surfaceUsed(1,:),surfaceUsed(2,:),surfaceUsed(3,:),'r.','MarkerSize',15);
hold off

% It probably only makes sense to interpolated values close to the original measurement locations!

