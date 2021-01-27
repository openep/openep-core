function [centroid, allCentroids] = tricentroid(tri)
% TRICENTROID calculates the centroid of a TriRep object.
% Usage:
%   [centroid, allCentroids] = tricentroid(tri)
% Where:
%   tri is a TriRep object
%   centroid is the centroid
%   allCentroids contains the centroid of each triangle in the triangulation
%
% TRICENTROID calculates the weighted average of all the centroids of all
% the triangles represented by tri.triangulation, with the weighting being
% the area of each triangle.
%
% Author: Nick Linton (2009)
% Modifications -
%       Steven Williams (2021): Modified to work with Triangulation class
%
% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------


% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

if isa(tri, 'TriRep')
    data.X = tri.X;
    data.Triangulation = tri.Triangulation;
elseif isa(tri, 'triangulation')
    data.X = tri.Points;
    data.Triangulation = tri.ConnectivityList;
else
    error('TRICENTROID: Input should be a TriRep object or a Triangulation object')
end

% first get the areas of the triangles
a = triarea(tri);

% now the centroids
allCentroids = 1/3 * (   data.X(data.Triangulation(:,1),:) ...
            + data.X(data.Triangulation(:,2),:) ...
            + data.X(data.Triangulation(:,3),:)          );
        
centroid = sum(a(:,ones(1,size(allCentroids,2))).*allCentroids,1) ./ sum(a,1);