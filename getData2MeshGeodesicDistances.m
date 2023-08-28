function [distances, isVertUsed, userdata] = getData2MeshGeodesicDistances(userdata, varargin)
% GETDATA2MESHGEODESICDISTANCES Calculates geodesic distances from data
% points to mesh points
%
% Usage:
%   distances = getData2meshGeodesicDistances(userdata)
%   [distances, isVertUsed] = getData2meshGeodesicDistances(userdata)
%   [distances, isVertUsed, userdata] = getData2meshGeodesicDistances(userdata)
%
% Where:
%   userdata  - an openEP data structure
%   distances  - an n x m array of distances, where n is the number of
%                electrode data points and m is the number of mesh nodes
%   isVertUsed - the vertices that are referenced by the triangulation
%
% GETDATA2MESHGEODESICDISTANCES repacks the original mesh to ensure that
% there are no vertices unreferenced by the triangulation. The vertices for
% which geodesic distances have been calculated are returned in isVertUsed.
% If three output argumetns are specified then a modified userdata
% structure containing the computed geodesic distances is also returned.
%
% See also: repack.m
%
% Author: Steven Williams (2016)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% int = openEpDataInterpolator
% distances = getData2MeshGeodesicDistances(userdata);
% d = int.interpolate(getVertices(userdata,'used',true), distances(1,:)', getVertices(userdata,'used',false));
% drawMap(userdata, 'type', 'bip', 'coloraxis', [0 120], 'data', d)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

global geodesic_library

if isfield(userdata.surface, 'geodesicDistances')
    distances = userdata.surface.geodesicDistances.distances;
    isVertUsed = userdata.surface.geodesicDistances.isVertUsed;
else
    % compute geodesic distances
    
    tic
    surfaceTriangulation = getMesh(userdata, 'type', 'triangulation');
    [surfaceTriangulation, isVertUsed] = repack(surfaceTriangulation);
    geodesic_library = 'geodesic_matlab_api';
    mesh = geodesic_new_mesh(surfaceTriangulation.Points, surfaceTriangulation.ConnectivityList);
    algorithm = geodesic_new_algorithm(mesh, 'exact');
    
    numPts = getNumPts(userdata);
    [~, startPts] = getElectrogramX(userdata, 'type', 'bip');
    
    distances = NaN(numPts, size(surfaceTriangulation.Points,1));
    for i = 1:numPts
        source_points{1} = geodesic_create_surface_point('vertex', i, startPts(i,:));
        geodesic_propagate(algorithm, source_points);
        [~, distances(i,:)] = geodesic_distance_and_source(algorithm);
    end
    
    geodesic_delete()
    
    T = toc;
    disp([ num2str(size(distances,2)) ...
          , ' geodesic distances from ' ...
          , num2str(size(distances,1)) ...
          , ' datapoints calculated in ' ...
          , num2str(T) ...
          , ' seconds.' ...
         ]);
     
     % store the data
     userdata.surface.geodesicDistances.distances = distances;
     userdata.surface.geodesicDistances.isVertUsed = isVertUsed;
     userdata.surface.geodesicDistances.computeTime = T;
end

clear geodesic_library

end