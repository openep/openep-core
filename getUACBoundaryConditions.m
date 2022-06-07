function [ivcIndices, svcIndices, csIndices, tvIndices,path, splittedMeshes] = getUACBoundaryConditions( userdata, rootMeshFile,varargin )
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
% figure
% hold on
% numberOfColorAllowed = 9;
% hSurf = drawMap(userdata, 'type', 'none');
% MV_Clipped = input('Is the Mitral Valve clipped?','s');

% if strcmp(MV_Clipped,'N')
%   [ newUserdata, newSurfaces ] = cutMesh( userdata )
%    newUserdata = setMesh(userdata, newSurfaces{1});
%    userdata = newUserdata;  
% end

figure
hSurf = drawMap(userdata, 'type', 'none');
% Get the free boundaries
[FF_temp, l, a, tr] = getAnatomicalStructures( userdata, 'plot', false);

for numFF = 1 : numel(FF_temp)
    FF_size(numFF) = numel(FF_temp{numFF});
end
[B sortedIndex] = sort(FF_size,'descend');
for i = 1 : 4 %numberOfColorAllowed 
    FF{i} = FF_temp{sortedIndex(i)};
end
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

%%%%%%%%%%%%%%%% This section commented by Ali ... needs to be cleaned later%%%%%%%%%%
% select a point near the SVC and a point near the IVC
% originalMesh = getMesh(userdata, 'type', 'triangulation');
% surfaceMesh = getMesh(userdata, 'type', 'triangulation', 'limittotriangulation', true);
% pp = PointPicker(originalMesh, hSurf, gca());
% disp('select a point near the superior vena cava and a point near the inferior vena cava');
% pause();
% 
% % get these points - the first should be the SVC
% pointIndices = pp.PointIndices;
% 
% % find the closest point in the SVC boundary for the SVC point
% 
% % So
% % pointIndices(1) is the SVC point, indexing into originalMesh.Points
% % pointIndices(2) is the IVC point, indexing into originalMesh.Points
% 
% % find the SVC point on the ring
% svcPointOnRing = findclosestvertex(originalMesh.Points(svcIndices,:), originalMesh.Points(pointIndices(1),:));
% svcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(svcIndices(svcPointOnRing),:));
%  
% ivcPointOnRing = findclosestvertex(originalMesh.Points(ivcIndices,:), originalMesh.Points(pointIndices(2),:));
% ivcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(ivcIndices(ivcPointOnRing),:));
% 
% plotTag(userdata, 'coord', surfaceMesh.Points(svcPointOnSurface,:), 'color', 'g');
% plotTag(userdata, 'coord', surfaceMesh.Points(ivcPointOnSurface,:), 'color', 'g');


% Initialise the geodesic library and algorithm

%%%%%%% Editted by Ali to able the code to make 3 paths
originalMesh = getMesh(userdata, 'type', 'triangulation');
surfaceMesh = getMesh(userdata, 'type', 'triangulation', 'limittotriangulation', true);

for i = 1 : 3
    
    if i == 1
        %select points near the SVC and near the IVC
        pp = PointPicker(originalMesh, hSurf, gca());
        disp('select a point near the superior vena cava and a point near the inferior vena cava');
        pause();
        pointIndices = pp.PointIndices;

        %find the SVC point on the ring
        svcPointOnRing = findclosestvertex(originalMesh.Points(svcIndices,:), originalMesh.Points(pointIndices(1),:));
        svcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(svcIndices(svcPointOnRing),:));
        %find the SVC point on the ring
        ivcPointOnRing = findclosestvertex(originalMesh.Points(ivcIndices,:), originalMesh.Points(pointIndices(2),:));
        ivcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(ivcIndices(ivcPointOnRing),:));
        % plot selected points on the mesh
        plotTag(userdata, 'coord', surfaceMesh.Points(svcPointOnSurface,:), 'color', 'g');
        plotTag(userdata, 'coord', surfaceMesh.Points(ivcPointOnSurface,:), 'color', 'g');
        % define the begining and the ending poits to calculate the geodesic path coneecting them
        startPoint = svcPointOnSurface;
        endpoint = ivcPointOnSurface;

    elseif i == 2
        %select points near the SVC and near the TV (tricuspid Valve)
        pp = PointPicker(originalMesh, hSurf, gca());
        disp('select a point near the superior vena cava and a point near the tricuspid valve');
        pause();
        pointIndices = pp.PointIndices;
        %find the SVC point on the ring
        svcPointOnRing = findclosestvertex(originalMesh.Points(svcIndices,:), originalMesh.Points(pointIndices(1),:));
        svcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(svcIndices(svcPointOnRing),:));
        %find the TV point on the ring
        TVPointOnRing = findclosestvertex(originalMesh.Points(tvIndices,:), originalMesh.Points(pointIndices(2),:));
        TVPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(tvIndices(TVPointOnRing),:));
        % plot selected points on the mesh
        plotTag(userdata, 'coord', surfaceMesh.Points(svcPointOnSurface,:), 'color', 'b');
        plotTag(userdata, 'coord', surfaceMesh.Points(TVPointOnSurface,:), 'color', 'b');
        % define the begining and the ending poits to calculate the geodesic path coneecting them
        startPoint = svcPointOnSurface;
        endpoint = TVPointOnSurface;
    else
        %select points near the IVC and near the TV (tricuspid Valve)

        pp = PointPicker(originalMesh, hSurf, gca());
        disp('select a point near the inferior vena cava and a point near the tricuspid valve');
        pause();
        pointIndices = pp.PointIndices;
        %find the IVC point on the ring
        ivcPointOnRing = findclosestvertex(originalMesh.Points(ivcIndices,:), originalMesh.Points(pointIndices(1),:));
        ivcPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(ivcIndices(ivcPointOnRing),:));
        %find the TV point on the ring
        TVPointOnRing = findclosestvertex(originalMesh.Points(tvIndices,:), originalMesh.Points(pointIndices(2),:));
        TVPointOnSurface = findclosestvertex(surfaceMesh, originalMesh.Points(tvIndices(TVPointOnRing),:));

        % plot selected points on the mesh
        plotTag(userdata, 'coord', surfaceMesh.Points(ivcPointOnSurface,:), 'color', 'k');
        plotTag(userdata, 'coord', surfaceMesh.Points(TVPointOnSurface,:), 'color', 'k');

        % define the begining and the ending poits to calculate the geodesic path coneecting them
        startPoint = ivcPointOnSurface;
        endpoint = TVPointOnSurface;
    end
    % calculate geodesic path
    initialiseGeodesic()
    mesh = geodesic_new_mesh(surfaceMesh.Points, surfaceMesh.ConnectivityList);
    algorithm = geodesic_new_algorithm(mesh, 'dijkstra');

    source_points{1} = geodesic_create_surface_point('vertex', startPoint, surfaceMesh.Points(startPoint,:));
    geodesic_propagate(algorithm, source_points);
    destination = geodesic_create_surface_point('vertex', endpoint, surfaceMesh.Points(endpoint,:));
    path{i} = geodesic_trace_back(algorithm, destination);
    [x,y,z] = extract_coordinates_from_path(path{i});
%     if i==3
%     path3_xyz = [x,y,z];
%     Idx = knnsearch(surfaceMesh.Points,path3_xyz,'dist', 'euclidean','k',10);
%     x = [x; surfaceMesh.Points(Idx(:),1)];
%     y = [y; surfaceMesh.Points(Idx(:),2)];
%     z = [z; surfaceMesh.Points(Idx(:),3)];
%     
%     end

plot3( x*1.001 ...
    ,y*1.001 ...
    ,z*1.001 ...
    ,'k-' ...
    , 'LineWidth', 2 ...
    , 'markersize', 10 ...
    );

geodesic_delete();
end

%%%%% Cut the mesh along three calculated paths
removeVerticeList = [];
for iPath = 1 : size(path,2)
    for iverticePath = 1 : size(path{iPath},1)
        removeVerticeList = [removeVerticeList;[path{iPath}{iverticePath}.x ...
            path{iPath}{iverticePath}.y  path{iPath}{iverticePath}.z]];
    end

end
fvcIn.vertices = userdata.surface.triRep.X;
fvcIn.faces = userdata.surface.triRep.Triangulation;
Idx = knnsearch(fvcIn.vertices,removeVerticeList,'k',10);
removeVerticeListThick = [removeVerticeList;fvcIn.vertices(Idx(:),:)];
fvcOut = removeVerticesPatch(fvcIn,removeVerticeList);
figure
trisurf(fvcOut.faces,fvcOut.vertices(:,1),fvcOut.vertices(:,2),fvcOut.vertices(:,3),'EdgeColor','none')
axis square tight equal vis3d
Clipping = 'off';
colormap('jet')
selpath = getLaplaceSolution (tvIndices, ivcIndices, svcIndices, path,userdata,rootMeshFile);
splittedMeshes = splitMesh(fvcOut.vertices,fvcOut.faces);
cd (sprintf('%s',selpath,'/OUTPUT_DIR'))

end

%end