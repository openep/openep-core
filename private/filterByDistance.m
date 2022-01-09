function newData = filterByDistance(data, loc, X, threshold)
% FILTERBYDISTANCE Is used to keep only data values within threshold
% distance of X.
%
% Usage:
%   newData = filterByDistance(data, loc, X, threshold)
% Where:
%   data        - the original data
%   loc         - locations of original data
%   X           - reference locations
%   threshold   - the distance threshold to apply
%
% FILTERBYDISTANCE works on linear distances
%
% Author: Steven Williams (2022)
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

% identify the vertices
% x1 are interpolated data i.e. all the vertices
% x0 real data i.e. co-ordinates of the mapping points
x1 = loc;
x0 = X;

id = knnsearch(x0, x1);
cPts = x0(id,:); %c for closest
d = distBetweenPoints(cPts, x1);
vtx = ones(size(d));
vtx(d>threshold) = 0;
vtx = logical(vtx);

newData = data;
newData(~vtx) = NaN;

end