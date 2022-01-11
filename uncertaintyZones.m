function [ points, weights ] = uncertaintyZones(userdata, loc, regionDefinition, varargin)
% UNCERTAINTYZONES identifies points of interest around locations of
% interest
%
% Usage:
%   [ points, weights ] = uncertaintyZones( userdata, loc, regionDefinition )
% Where:
%   userdata  - see importcarto_mem
%   loc  - the co-ordinates of the points we are interested in, n x 3
%   regionDefinition - structure with the following fields:
%       .shape - the shape, only sphere is implemented
%       .params - parameters describing the shape, only radius is applicable
%   points - logical array of size n x m
%   weights - numerical array of size n x m
%
% UNCERTAINTYZONES accepts the following parameter-value pairs
%   'type'              {'surface'}|'electrogram'
%   'weightalgorithm'   {'linear'}
%
% UNCERTAINTYZONES is used to identify points of interest around locations
% of interest. For example, UNCERTAINTYZONES can be used to access data
% such as electrogram voltage data within a region around each mapping
% point.
%
% The function works by calculating a radius aroudn each location of
% interest (loc), and identifying every data point that falls within that
% sphere. A weighting is calculated for each data point of interest based
% on the weighting algorithm, currently the only option is _linear_.
%
% UNCERTAINTYZONES can work on surface data or electrogram data. In the
% case of surface data, the mesh nodes are taken as the set of all possible
% data locations. In the case of electrogram data, the electrogram
% recordings locations are taken as the set of all possible data locations.
%
% The results are returned as two arrays, points and weights. The size of
% both is n x m, where n is the number of locations queried, and m is the
% number of all data locations available. Points takes values of 0 at all
% data locations outwith the scope of the region definition and 1 at all
% locations within the scope of the region definition. Weights has a
% numberical value representing the weights at all locations, in the range
% [0 1]. By definition, any data locations outwith the region definition
% will have a weight of zero.
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

% Parse inputs
nStandardArgs = 3; % UPDATE VALUE
type = 'map';
weightAlgorithm = 'linear';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'type'
                type = varargin{i+1};
            case 'weightalgorithm'
                weightAlgorithm = varargin{i+1};
        end
    end
end

% Create default regionDefinition
if isempty(regionDefinition)
    regionDefinition.shape = 'sphere';
    regionDefinition.params = 5;
    disp('Default region definition used');
end

% get the potential data positions
switch lower(type)
    case 'map'
        X = getVertices(userdata); % todo - how do we handle points not referenced by the triangulation
    case 'electrogram'
        X = getElectrogramX(userdata);
end

% identify the points of interest
m = size(X,1);
n = size(loc,1);
points = zeros(n, m);
weights = zeros(n,m);

switch lower(regionDefinition.shape)
    case 'sphere'
        radius = regionDefinition.params(1);
        for n = 1:size(loc,1)
            distances = distBetweenPoints(loc(n,:),X,'method','linear');
            points(n,distances<radius) = 1;

            temp = distances/radius-1;
            temp(temp>0) = 0;
            weights(n,:) = abs(temp);
        end
end

points = points';
weights = weights';

end







