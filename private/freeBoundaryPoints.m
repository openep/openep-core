function [coords] = freeBoundaryPoints(FF, tr)
% FREEBOUNDARYPOINTS returns the co-ordinates of the vertices on the free
% boundary, FF, of triRep, tr
%
% Usage:
%   [ coords ] = freeBoundaryPoints( FF, tr )
% Where:
%   tr  - a triRep object
%   FF  - the free boundary facets, see freeBoundary.m
%   coords  - the co-ordinates of the vertices on the free boundary FF, of
%             the form:
%         [ x_1  y_1  z_1 ]
%         [ x_2  y_2  z_2 ]
%         [ ...  ...  ... ]
%         [ x_n  y_n  z_n ]
%
% FREEBOUNDARYPOINTS Detailed description goes here
%
% See also: FREEBOUNDARIES, DRAWFREEBOUNDARIES
%
% Author: Steven Williams (2020)
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

coords = tr.X(FF(:,1),:);

end

