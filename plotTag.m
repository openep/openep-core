function h = plotTag( userdata, varargin )
% PLOTTAG Plots tag(s) on the current map
%
% Usage:
%   h = plotTag( userdata, varargin )
% Where:
%   userdata  - see importcarto_mem
%   h(i) - is an array of handles referencing the plotted surfaces
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
        end
    end
end
% TODO add error checking for input parsing

hold on
if ~isempty(coord)
    % plot coordinates
    for i = 1:size(coord,1)
        h(i) = plotsphere(coord(i,1), coord(i,2), coord(i,3), color, r, 16);
    end
    
end

if ~isempty(pointnum)
    % plot pointnums
    coord = userdata.electric.egmX(pointnum,:);
    for i = 1:numel(pointnum)
        h(i) = plotsphere(coord(i,1), coord(i,2), coord(i,3), color, r, 16);
    end
end

set(h, 'edgecolor', 'none');

end
