function [ vtx ] = getVerticesNearMappingPoints( userdata, distanceThreshold, varargin )
% GETVERTICESNEARMAPPINGPOINTS identifies mesh vertices which are within
% distanceThreshold of an actual mapping point
%
% Usage:
%   [ vol ] = getVerticesNearMappingPoints( userdata, distanceThreshold )
% Where:
%   userdata  - see importcarto_mem
%   distanceThreshold  - the distance threshold to be applied
%   vtx  - logical array of size length(userdata.surface.triRep.X) which
%          identifies the vertices within distanceThreshold of an actual
%          mapping point.
%
% GETVERTICESNEARMAPPINGPOINTS accepts the following parameter-value pairs
%   'method'     {'linear'}|'geodesic'    TODO: N.B. geodesic not yet implemented
%
% GETVERTICESNEARMAPPINGPOINTS Can be used to identify mesh vertices which
% are within a particular distanceThreshold of an original mapping point.
% For example, the function can identify all the mesh vertices within 10mm
% of a local activation time mapping point. This is useful in order to
% filter interpolated data by only those points that are close to points of
% actual data. Two methods are provided for calculating distances:
% Euclidean distances and geodesic distances. Note that geodesic distances
% is not yet implemented.
%
% Author: Steven Williams (2021)
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

nStandardArgs = 2; % UPDATE VALUE
method = 'linear';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
        end
    end
end

% identify the vertices
% x1 are interpolated data i.e. all the vertices
% x0 real data i.e. co-ordinates of the mapping points
x1 = getVertices(userdata);
x0 = userdata.electric.egmSurfX;

if strcmpi(method, 'geodesic')
    error('OPENEP/GETVERTICESNEARMAPPINGPOINTS: Geodesic distances not yet implemented');
    %TODO: Geodesic functionality needs to be added to distBetweenPoints.m; when working, this function will work too
end

id = knnsearch(x0, x1);
cPts = x0(id,:); %c for closest
d = distBetweenPoints(cPts, x1, 'method', method);
vtx = ones(size(d));
vtx(d>distanceThreshold) = 0;
vtx = logical(vtx);

end