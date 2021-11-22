function [D, pathCoordinates] = distBetweenPoints(A, B, varargin)
% DISTBETWEENPOINTS Returns the distance from A to B.
% Usage:
%   D = DISTBETWEENPOINTS(A, B)
% Where:
%   A - is the first point(s)
%   B - is the second point(s)
%   D - distance between the points
%   pathCoordinates - co-ordinates of the geodesic path, empty if 'method'
%                     is 'linear'
%
% DISTBETWEENPOINTS accepts the folloiwng parameter-value pairs
%   'method'    {'linear'} | 'geodesic'
%       - Specifies whether to calcualte linear or geodesic distances
%   'userdata'  {[]} | userdata
%       - Specifies the OpenEP userdata structure for geodesic calculations
%
% DISTBETWEENPOINTS returns the distance from A to B. A and B are specified
% as row vectors [x, y, z] or matrices, with rows representing different
% points. If npoints in A and B are different A must specify one and only 
% one point.
%
% Author: Steven Williams (2013)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% set up global variables
global geodesic_library;
geodesic_library = 'geodesic_matlab_api';

% parse input arguments
nStandardArgs = 3;
method = 'linear';
warn = false;
userdata = [];
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'verbose'
                warn = varargin{i+1};
            case 'userdata'
                userdata = varargin{i+1};
        end
    end
end

if size(A,2) ~= size(B,2)
    error('DISTBETWEENPOINTS: dimensions of points in A and B must be equal');
end
if size(A,1) ~= size(B,1)
    if size(A,1) ~= 1
        error('DISTBETWEENPOINTS: if npoints in A and B are different A must specify one and only one point');
    else
        A = repmat(A, size(B,1), 1);
        if warn
            warning('DISTBETWEENPOINTS: A has been replicated size(B,1) times');
        end
    end
end

switch lower(method)
    case 'linear'
        diffsq = (A - B).^2;
        D = sqrt(sum(diffsq, 2));
        pathCoordinates = [];

    case 'geodesic'
        if isempty(userdata)
            error('OPENEP/DISTBETWEENPOINTS: you must specify userdata for geodesic distance calculations')
        end
        % calculate the geodesic distance (repacking and reducepatch
        % necessary to prevent ?memory problem in exact geodesic)
        V = userdata.surface.triRep.X;
        F = userdata.surface.triRep.Triangulation;
        [faces,vertices] = reducepatch(F,V,size(F,1));
        
        % now find the closest vertex index to A and B
        P1new = findclosestvertex(vertices, A);
        P2new = findclosestvertex(vertices, B);
        
        % geodesic algorithm
        mesh = geodesic_new_mesh(vertices,faces);
        algorithm = geodesic_new_algorithm(mesh, 'exact');
        source_point = {geodesic_create_surface_point('vertex',P1new,vertices(P1new,:))};
        geodesic_propagate(algorithm, source_point); % the most time-consuming step
        
        % find a shortest path from source to target
        destination = geodesic_create_surface_point('vertex',P2new,vertices(P2new,:));
        path = geodesic_trace_back(algorithm, destination);
        
        % find distances
        [x,y,z] = extract_coordinates_from_path(path);
        D = sum(sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2));

        pathCoordinates = [x y z];
        
        %delete all meshes and algorithms
        geodesic_delete;

end
