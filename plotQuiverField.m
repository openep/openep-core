function h = plotQuiverField(userdata, varargin)
% PLOTQUIVERFIELD Plots vector data on a surface.
%
% Usage:
%   h = plotQuiverField(userdata)
%   h = plotQuiverField(userdata, cvX, u)
%
% Where:
%   userdata - an OpenEP data structure
%   h        -  a quivergroup handle
%
% PLOTQUIVERFIELD accepts the following parameter-value pairs
%   'rbfoptions'     see doCvMapping_RadialBasis
%   'X'              the starting positions, size n x 3
%   'u'              the vectors to be plotted, size n x 3
%
% PLOTQUIVERFIELD plots vector data on a surface. For example
% PLOTQUIVERFIELD can be used to plot activation directions on a local
% activation map, which is the default behaviour. Optionally, vectors can 
% be specified as two matrices, X, the starting locations and u the vectors
% to be plotted. This can be used for example, to display other forms of
% vector data as a quiver field in OpenEP.
%
% Author: Steven Williams (2021)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% % Example 1: create a quiver plot using default settings
% load openep_dataset_1
% plotQuiverField(userdata);
%
% % Example 2: specify the vectors to be plotted
% load openep_dataset_1
% [~, cvX, ~, u] = getConductionVelocity(userdata, 'method', 'gradient');
% plotQuiverField(userdata, 'X', cvX, 'u', u);
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse input arguments
nStandardArgs = 1; % UPDATE VALUE
X = [];
u = [];
xSpecified = false;
uSpecified = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'x'
                X = varargin{i+1};
                xSpecified = true;
            case 'u'
                u = varargin{i+1};
                uSpecified = true;
        end
    end
end

if xSpecified & ~uSpecified
    error('OPENEP/PLOTQUIVERFIELD: If X is specified, then u must also be specified');
end
if ~xSpecified & uSpecified
    error('OPENEP/PLOTQUIVERFIELD: If u is specified, then X must also be specified');
end

% calculate the conduction vectors, if necessary
if ~xSpecified && ~uSpecified
    [~, X, ~, u] = getConductionVelocity(userdata, 'method', 'radialbasis');
end

% get the surface normals
[normals, userdata] = getNormals(userdata);
iV = findclosestvertex(getMesh(userdata, 'trirep'), X);
normals = normals(iV,:);

% project the velocities in tangent direction
projU=projectVertexVectors(X,u,normals);

% compute unit-vectors:
uv=createUnitVectors(projU);

% plot
hold on
h = quiver3(X(:,1),X(:,2),X(:,3),uv(:,1),uv(:,2),uv(:,3) ...
    , 'color', 'k' ...
    , 'linewidth', 1 ...
    );

end