function h = drawFreeBoundaries(FF, tr, varargin)
% DRAWFREEBOUNDARIES Draws the boundaries at the edge of a trirep
%
% Usage:
%   h = drawFreeBoundaries(FF, tr)
%
% Where:
%   FF  - a cell array of connected free boundary facets, see freeBoundary.m
%   tr  - is a triRep object
%   h   - a cell array of shandle to the line
%
% DRAWFREEBOUNDARIES draws the free boundaries around a trirep.
%
% Author: Steven Williams (2020)
% Modifications -
%
% See also: FREEBOUNDARIES, FREEBOUNDARYPOINTS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

for i = 1:numel(FF)
    h{i} = local_draw(FF{i}, tr, colorBrewer(i));
end

    function h = local_draw(rimEdges, tr, col)
        n = size(rimEdges,1);
        x = nan(3*n,1);
        y = nan(size(x));
        z = nan(size(x));
        
        for j = 0:(n-1)
            loc1 = tr.X(rimEdges(j+1,1),:);
            loc2 = tr.X(rimEdges(j+1,2),:);
            
            x(3*j+1) = loc1(1);
            y(3*j+1) = loc1(2);
            z(3*j+1) = loc1(3);
            
            x(3*j+2) = loc2(1);
            y(3*j+2) = loc2(2);
            z(3*j+2) = loc2(3);
        end
        
        h = line();
        set(h ...
            , 'XData', x ...
            , 'YData', y ...
            , 'ZData', z ...
            , 'LineWidth', 2 ...
            , 'Color', col ...
            , 'visible', 'on' ...
            );
    end

end