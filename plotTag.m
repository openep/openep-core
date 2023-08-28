function [h vi] = plotTag( userdata, varargin )
% PLOTTAG Plots tag(s) on the current map
%
% Usage:
%   h = plotTag( userdata, varargin )
%   
% Where:
%   userdata  - see importcarto_mem
%   h - an array of handles representing the plotted surface
%   vi - indices of the surface projection if 'type' is 'surf' or 'both'
%
% PLOTTAG accepts the following parameter-value pairs
%   'coord'     {[x y x]}
%       - A set of x,y,z coords where size(coords) = nx3 where n is the
%       number of tags to plot
%   'pointnum'  [p1, p2, ... pn]
%       - An array of size nx1 where n is the number of tags to plot
%   'color'     {'r'}|'g'|'b'|'p'|'o'|'y'
%       - The color of the tag to draw
%   'size'
%       - The size of the tag to draw
%   'type'      {'3d'}|'surf'|'both'
%
% PLOTTAG Plots tag(s) on the current map
%
% TODO: Complete implementation of pointTag.m
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

nStandardArgs = 1; % UPDATE VALUE
coord = [];
pointnum = [];
color = 'r';
r = 3;
type = '3d';
if nargin > nStandardArgs
    for i = 1:2:nargin-1
        switch varargin{i}
            case 'coord'
                coord = varargin{i+1};
            case 'pointnum'
                pointnum = varargin{i+1};
            case 'color'
                color = varargin{i+1};
            case 'size'
                r = varargin{i+1};
            case 'type'
                type = varargin{i+1};
        end
    end
end
% TODO add error checking for input parsing

% set up the constants
SURFDIAMETER = r + 0.5;
SURFTHICKNESS = 1.5;
NUMELEM = 16;

hold on
if ~isempty(pointnum)
    coord = [coord; userdata.electric.egmX(pointnum,:)];
end

vi = [];
if strcmpi(type, 'surf') || strcmpi(type, 'both')
    if ~isfield(userdata.surface, 'normals')
        N = getNormals(userdata);
    else
        N = userdata.surface.normals;
    end
    if ~isempty(pointnum)
        surfCoord = userdata.electric.egmSurfX(pointnum,:);        
        vi = findclosestvertex(getMesh(userdata, 'type', 'triangulation', 'limittotriangulation', false), surfCoord);
    else
        mesh = getMesh(userdata, 'type', 'triangulation', 'limittotriangulation', false);
        vi = findclosestvertex(mesh, coord);
        surfCoord = mesh.Points(vi,:);
    end
    endPoints = surfCoord + (N(vi,:) * SURFTHICKNESS);
    
end

h = [];
for i = 1:size(coord,1)
    
    if strcmpi(type, '3d') || strcmpi(type, 'both')
        h(end+1) = plotsphere(coord(i,1), coord(i,2), coord(i,3), color, r, NUMELEM);
    end
    if strcmpi(type, 'surf') || strcmpi(type, 'both')
        ind = length(h);
        [h(ind+1) h(ind+2) h(ind+3)] = drawCylinder([surfCoord(i,1), surfCoord(i,2), surfCoord(i,3) ...
            , endPoints(i,1), endPoints(i,2), endPoints(i,3) ...
            , SURFDIAMETER], NUMELEM ...
            , 'FaceColor', color ...
            );
        set(h(2:3), 'DiffuseStrength', 0.5);
    end
    
end

if strcmpi(type, '3d') || strcmpi(type, 'both')
    set(h, 'edgecolor', 'none');
end

end
