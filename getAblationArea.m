function [ablArea, isAblated, trAbl] = getAblationArea(userdata, varargin)
% GETABLATIONAREA Calculates the area of a chamber which has been ablated
%
% Usage:
%   [ablArea, isAblated, trAbl] = getAblationArea( userdata )
% Where:
%   userdata - see importcarto_mem.m
%   ablArea - the total area of the chamber that has been ablated
%   isAblated - indexes into getMesh(userdata).Triangulation and
%               indicates whether a particular triangle is considered
%               ablated (1) or not (0).
%   trAbl - a Triangulation of the ablated tissue
%
% GETABLATIONAREA accepts the following parameter-value pairs
%   'method'     {'tags'}|'grid'
%       - specifies whether to calculate area based on the ablation tags or
%         the ablation grid
%   'radius'    {5}|double
%       - specifies the radius around each ablation tag to consider ablated
%   'thresholdmethod'    {'on'}|'off'
%       - NOT YET IMPLEMENTED
%   'thresholdvalue'
%       - NOT YET INMPLEMENTED
%
% GETABLATIONAREA Requires a userdata structure which contains .rfindex as
% its input, which can be created using importvisitag.m
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
%
% See also, importvisitag.m, plotAblationArea.m, plotVisitags.m

nStandardArgs = 1; % UPDATE VALUE
method = 'tags';
r = 5;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                toplot = varargin{i+1};
            case 'radius'
                r = varargin{i+1};
        end
    end
end
%TODO add error checking for the input param/value pairs

switch method
    case 'tags'
        X = userdata.rfindex.tag.X;
        
        % find the closest point on the surface to each tag
        iV = findclosestvertex(getMesh(userdata), X);
        vertices = getMesh(userdata).X(iV,:);
        
        % find the centroids of all the triangles in the triangulation
        [~,allCentroids] = tricentroid(getMesh(userdata));
        
        % find all the triangles whose centres are within r radius of
        isAblated = false(size(getMesh(userdata).Triangulation,1),1);
        for i = 1:length(vertices)
            distances = distBetweenPoints(vertices(i,:), allCentroids);
            isAblatedByThisLesion = distances<r;
            isAblated(isAblatedByThisLesion) = true;
        end
        trAbl = triangulation(getMesh(userdata).Triangulation(isAblated,:) ...
            , getMesh(userdata).X(:,1) ...
            , getMesh(userdata).X(:,2) ...
            , getMesh(userdata).X(:,3) ...
            );
        ablArea = sum(triarea(trAbl))/100;
        
    case 'grid'
        warning('OPENEP/GETABLATIONAREA: Grid method not yet implemented')
        
end

end