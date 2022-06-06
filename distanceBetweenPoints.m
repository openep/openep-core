function distance = distanceBetweenPoints(userdata, P1, P2, varargin)
% DISTANCEBETWEENPOINTS Returns the distance from A to B.
% Usage:
%   distance = distanceBetweenPoints(userdata, A, B)
% Where:
%   userdata    - see importcarto_mem
%   P1          - is the first point
%   P2          - is the second point
%
% DISTBETWEENPOINTS accepts the following parameter-value pairs
%   'method'    {'linear'} | 'geodesic'
%       - Specifies whether to calcualte linear or geodesic distances
%   'plot'      {false} | true
%       - Specifies whether to draw a figure
%
% DISTBETWEENPOINTS returns the distance from A to B. A and B are specified
% as row vectors [x, y, z] or matrices, with rows representing different
% points. If npoints in A and B are different A must specify one and only
% one point.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
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

% parse input arguments
nStandardArgs = 3;
method = 'linear';
plot = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'plot'
                plot = varargin{i+1};
        end
    end
end

switch method
    case 'linear'
        % get the co-ordinates of the points
        A = userdata.electric.egmX(P1,:);
        B = userdata.electric.egmX(P2,:);
        
        % calculate the linear distance
        diffsq = (A - B).^2;
        distance = sqrt(sum(diffsq, 2));
        
        % plot a figure
        if plot
            hSurf = drawMap(userdata, 'type', 'none', 'orientation', 'pa');
            set(hSurf, 'facealpha', .5);
            hold on
            plotTag(userdata, 'coord', A);
            plotTag(userdata, 'coord', B);
            line([A(1) B(1)], [A(2) B(2)], [A(3) B(3)], 'linewidth', 4);
        end
        
    case 'geodesic'
        % get the co-ordinates of the surface points
        A = userdata.electric.egmSurfX(P1,:);
        B = userdata.electric.egmSurfX(P2,:);
        
        % calculate the geodesic distance (repacking to remove points
        % not referenced in the triangulation)
        tr = getMesh(userdata, 'limitToTriangulation', true);
        vertices = tr.X;
        faces = tr.Triangulation;
        
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
        distance = sum(sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2));
        
        %delete all meshes and algorithms
        geodesic_delete;
        
        if plot
            hSurf = drawMap(userdata, 'type', 'none', 'orientation', 'pa');
            set(hSurf, 'facealpha', .5);
            hold on
            plotTag(userdata, 'coord', A);
            plotTag(userdata, 'coord', B);
            plot3(x*1.001,y*1.001,z*1.001,'k-','LineWidth',2);    %plot path
        end
end

end
