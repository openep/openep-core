function distance = distanceBetweenPoints(userdata, P1, P2, varargin)
% DISTANCEBETWEENPOINTS Returns the distance between electrogram recording 
% locations P1 and P2.
%
% Usage:
%   distance = distanceBetweenPoints(userdata, A, B)
% Where:
%   userdata    - see importcarto_mem
%   P1           - is the index of the first electrogram point(s)
%   P2           - is the index of the second electrogram point(s)
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

% set up global variables
global geodesic_library;
geodesic_library = 'geodesic_matlab_api';

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

        distance = distBetweenPoints(A, B, 'method', 'linear');
        
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

        [distance, pathCoordinates] = distBetweenPoints(A, B, 'method', 'geodesic', 'userdata', userdata);
        
        if plot
            hSurf = drawMap(userdata, 'type', 'none', 'orientation', 'pa');
            set(hSurf, 'facealpha', .5);
            hold on
            plotTag(userdata, 'coord', A);
            plotTag(userdata, 'coord', B);
            x = pathCoordinates(:,1);
            y = pathCoordinates(:,2);
            z = pathCoordinates(:,3);
            plot3(x*1.001,y*1.001,z*1.001,'k-','LineWidth',2);    %plot path
        end
end

end
