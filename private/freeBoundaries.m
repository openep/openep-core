function [FF, l] = freeBoundaries(tr)
% FREEBOUNDARIES returns a cell array of connected free boundaries
%
% Usage:
%   [ FF, l ] = freeBoundaries( tr )
% Where:
%   tr  - a triRep object
%   FF  - a cell array of connected free boundary facets, see freeBoundary.m
%   l  - the lengths of the free boundaries in FF
%
% FREEBOUNDARIES Detailed description goes here
%
% Author: Steven Williams (2020)
% Modifications -
%
% See also: DRAWFREEBOUNDARIES, FREEBOUNDARYPOINTS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% Find all the free boundarie facets
FFall = freeBoundary(tr);

% if the edges are continuous then FF(1,i+1) == FF(2,i), for
% 1<i<length(FF). So, we first offset the columns in FF and check this
% equality
testset = NaN(size(FFall));
testset(:,1) = FFall(:,1);
testset(:,2) = [FFall(end,2); FFall(1:end-1,2)];
iStart = find(diff(testset,1,2));

% iStart are the indices of the start facets of the free boundaries. The
% number of connected free boundaries (i.e. the number of holes in the
% mesh) is given by numel(iStart)
if isempty(iStart)
    FF{1} = FFall;
    coords = freeBoundaryPoints(FF{1}, tr);
    l = lineLength(coords);
else
    for i = 1:numel(iStart)
        if i<numel(iStart)
            FF{i} = FFall(iStart(i):iStart(i+1)-1,:);
            coords = freeBoundaryPoints(FF{i},tr);
            l(i) = lineLength(coords);
        else
            FF{i} = FFall(iStart(i):end,:);
            coords = freeBoundaryPoints(FF{i},tr);
            l(i) = lineLength(coords);
        end
    end
end       


end

