function h = drawFreeBoundary(tr, col, varargin)
% DRAWFREEBOUNDARY Draws a boundary at the edge of a trirep.
% Usage:
%   h = drawFreeBoundary(tr, col)
%   h = drawFreeBoundary(tr, col, false)
% Where:
%   tr is the triangulation
%   col is the desired color in rgb matrix
%   varargin can be true or false (to return xyz and draw nothing)
%   h is the handle to the line (or a structure h.X, h.Y, h.Z)
%
% DRAWFREEBOUNDARY draws the free boundary around a trirep. Passing false
% as the third argument means that the boundary is not actually drawn but
% x, y and z are returned in the structure h
%
% Author: Steven Williams (2014)
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

drawline = true;
if nargin==3
    drawline = varargin{1};
end

rimEdges= freeBoundary(tr);

n = size(rimEdges,1);
x = nan(3*n,1);
y = nan(size(x));
z = nan(size(x));

if isa(tr, 'TriRep')
    X = tr.X;
elseif isa(tr, 'triangulation')
    X = tr.Points;
end

for i = 0:(n-1)
    loc1 = X(rimEdges(i+1,1),:);
    loc2 = X(rimEdges(i+1,2),:);
    
    x(3*i+1) = loc1(1);
    y(3*i+1) = loc1(2);
    z(3*i+1) = loc1(3);
    
    x(3*i+2) = loc2(1);
    y(3*i+2) = loc2(2);
    z(3*i+2) = loc2(3);
end

if drawline
h = line();
set(h ...
    , 'XData', x ...
    , 'YData', y ...
    , 'ZData', z ...
    , 'LineWidth', 2 ...
    , 'Color', col ...
    , 'visible', 'on' ...
    );
else
    h.X = x;
    h.Y = y;
    h.Z = z;
end

end