function C = getCentreOfMass( userdata, varargin )
% GETCENTREOFMASS Returns the centre of mass of the anatomical model 
% defined in userdata
%
% Usage:
%   C = getCentreOfMass( userdata, varargin )
% Where:
%   userdata   - see importcarto_mem
%   C - the Cartesian co-ordinates of the centre of mass
%
% GETCENTREOFMASS accepts the following parameter-value pairs
%   'plot'     {false}|true
%
% GETCENTREOFMASS calculates the centre of mass of the userdata by
% accessing a closed surface via the OpenEP function getClosedSurface.m
% before using centroidOfPolyhedron.m to calculate the centre of mass. The
% function centroidOfPolyhedron.m was written by Isfandiyar RASHIDZADE,
% available through the Mathworks FileExchange:
%     https://www.mathworks.com/matlabcentral/fileexchange/63614-centroid-of-triangulated-polyhedron
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% See also FREEBOUNDARYPOINTS, GETANATOMICALSTRUCTURES
%
% Info on Code Testing:
% ---------------------------------------------------------------
% C = getCentreOfMass( userdata, 'plot', true )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse input
nStandardArgs = 1; % UPDATE VALUE
plot = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'plot'
                plot = varargin{i+1};
        end
    end
end

% get a trirep
trMesh = getClosedSurface(userdata);

% get all the co-ordinates
C = centroidOfPolyhedron(trMesh.X, trMesh.Triangulation);

% plot
if plot
    figure
    hS = trisurf(trMesh);
    set(hS ...
        , 'facecolor', [.5 .5 .5] ...
        , 'edgecolor', 'none' ...
        , 'facealpha', .2 ...
        );
    set(gcf, 'color', 'w' ...
        );
    axis off equal vis3d;
    hold on;
    plotsphere(C(1), C(2), C(3), 'r', 5, 16);
end

end