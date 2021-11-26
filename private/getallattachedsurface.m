function [result, goodTriangles, goodPoints] = getallattachedsurface(inputSurface, inputStart, startType, varargin)
% GETALLATTACHEDSURFACE Returns all of the surface continuous with start.
% Usage:
%   newTriRep = getallattachedsurface(oldTriRep, start, starttype)
%   newFaceIndices  = getallattachedsurface(oldFaces,  start, starttype)
%   [ ... , newFaceIndices] = getallattachedsurface( ... )
% Where:
%   newTriRep.X contains only data points that are continuous with start
%            .Triangulation is the corresponding triangulation
%   newFaceIndices -
%           oldFaces(newFaceIndices,:)  or  newTriRep(newFaceIndices,:)
%           gives the triangulation that only refers to data points
%           continuous with start.
%   start is either an index representing a vertex or one of the faces -
%   startType is either 'Face' or 'Vertex'
%
% GETALLATTACHEDSURFACE accepts the following parameter-value pairs
%   'verbose'     {true}|false
%
% Author: Nick Linton (2010)
% Modifications -
%   Steven Williams (2021) - added verbose option
%
% Info on Code Testing:
% ---------------------
% test code
% ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 3; % UPDATE VALUE
verbose = true;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'verbose'
                verbose = varargin{i+1};
        end
    end
end

if isa(inputSurface, 'triangulation')
    tri = inputSurface;
elseif isa(inputSurface, 'TriRep')
    tri = triangulation(inputSurface.Triangulation, inputSurface.Points);
else
    if size(inputSurface,2) ~= 3
        error('GETALLATTACHEDSURFACE: oldFaces should be a n*3 matrix')
    end
    tri = triangulation(inputSurface, zeros(max(inputSurface(:)),2));  %create a TriRep with dummy X data - this will let us use the built in functions for TriRep
end

switch lower(startType)
    case 'face'
        goodTriangles = inputStart;
    case 'vertex'
        goodTriangles = vertexAttachments(tri,inputStart);
        goodTriangles = goodTriangles{1};
        goodTriangles = goodTriangles(:);
    otherwise
        error('GETALLATTACHEDSURFACE: incorrect string for "startType"')
end

% Each face has 3 edges. There will be <= 3 neighbouring triangles for each
% face.
if verbose
    warning('if a triangle edge has 3 face attachments then this will not work')
end
ngbrs = neighbors(tri);
ngbrs = [ (1:size(ngbrs,1))' , ngbrs ]; %I'm including the index to the original triangles as well

% set the variables as large arrays to try and reserve contiguous memory
% for them (not sure that this approach works).
temp = goodTriangles;
goodTriangles = zeros(size(tri,1),1);
oldResult = nan(size(goodTriangles));
goodTriangles = temp;

while ~isequal(oldResult,goodTriangles)
    oldResult = goodTriangles;
    newGoodTriangles = ngbrs(goodTriangles,:);
    goodTriangles = unique(newGoodTriangles(isfinite(newGoodTriangles)));
end

if isa(inputSurface, 'triangulation')
    newTriangulation = tri.ConnectivityList(goodTriangles,:);
    [goodPoints, ~, newTriangulation(:)] = unique(newTriangulation(:));
    newX = tri.Points(goodPoints,:);
    result = triangulation(newTriangulation, newX);
elseif isa(inputSurface, 'TriRep')
    newTriangulation = tri.ConnectivityList(goodTriangles,:);
    [goodPoints, ~, newTriangulation(:)] = unique(newTriangulation(:));
    newX = tri.Points(goodPoints,:);
    result = TriRep(newTriangulation, newX); %#ok<DTRIREP>
else
    result = tri.ConnectivityList(goodTriangles,:);
    goodPoints = unique(result(:));
end
