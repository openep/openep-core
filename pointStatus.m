function [inoutpts, meshpts] = pointStatus( userdata, varargin )
% POINTSTATUS Returns the status of points relevant to userdata
%
% Usage:
%   [inoutpts, meshpts] = pointStatus( userdata )
% Where:
%   userdata  - see importcarto_mem
%   inoutpts  - whether points are internal (logical(1)) or external
%               (logical(0)) to the triangulation in userdata
%   meshpts   - whether points in the triangulation in userdata are
%               referenced in the triangulation (logical(1)) or not
%               (logical(0))
%
% POINTSTATUS accepts the following parameter-value pairs:
%   'tol' 0.1 | double
%       - The distance threshold within which points are considered to be
%         internal or external to the triangulation
%   'plot'     {false}|true
%       - Specify whether to plot the results
%
% POINTSTATUS depends on the package INPOLYHEDRON. See:
%   https://uk.mathworks.com/matlabcentral/fileexchange/37856-inpolyhedron-are-points-inside-a-triangulated-volume
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% load('/Users/Steven/Desktop/pepr_working_dir/503.mat')
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 1; % UPDATE VALUE
plot = false;
tol = 0.1;
if nargin > nStandardArgs
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'plot'
                plot = varargin{i+1};
            case 'tol'
                tol = varargin{i+1};
        end
    end
end

% First deal with mapping points
FV.vertices = getVertices(userdata);
FV.faces = getFaces(userdata);
QPTS = userdata.electric.egmX;
inoutpts = inpolyhedron(FV, QPTS, 'tol', tol);

% Now deal with surface points
[pts, meshpts] = getVertices(userdata);

if plot
    drawMap(userdata, 'type', 'none')
    hold on;
    
    % plot the inoutpts
    plot3(QPTS(inoutpts,1), QPTS(inoutpts,2), QPTS(inoutpts,3), '.g', 'markersize', 10);
    plot3(QPTS(~inoutpts,1), QPTS(~inoutpts,2), QPTS(~inoutpts,3), '.r', 'markersize', 10);
    
    % plot the meshpts
    plot3(pts(meshpts,1), pts(meshpts,2), pts(meshpts,3), '.k');
    plot3(pts(~meshpts,1), pts(~meshpts,2), pts(~meshpts,3), '.', 'color', [.5 .5 .5], 'markersize', .5);
end

end