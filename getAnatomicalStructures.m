function [FF, l, a, tr] = getAnatomicalStructures( userdata, varargin )
% GETANATOMICALSTRUCTURES Returns the free boundaries (anatomical 
% structures) described in userdata
%
% Usage:
%   [FF, l, a, tr] = getAnatomicalStructures( userdata, varargin )
% Where:
%   userdata   - see importcarto_mem
%   FF         - see TriRep/freeBoundary, cell array
%   l          - array of lengths (perimeters) of each anatomical structure
%   a          - an array of areas of each anatomical structure
%   tr         - cell array of triangulations of each anatomical structure
%
% GETANATOMICALSTRUCTURES accepts the following parameter-value pairs
%   'plot'     {false}|true
%
% GETANATOMICALSTRUCTURES identifies all the anatomical structures of a 
% given data set. Anatomical structures are boundary regions that have been 
% added to an anatomical model in the clinical mapping system. For example, 
% with respect of left atrial ablation, anatomical structures may represent 
% the pulmonary vein ostia, mitral valve annulus or left atrial appendage 
% ostium.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% See also FREEBOUNDARYPOINTS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse input
nStandardArgs = 1; % UPDATE VALUE
plot = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'plot'
                plot = varargin{i+1};
        end
    end
end

% get a trirep
trMesh = getMesh(userdata);

% get all the free boundaries
[ FF, l ] = freeBoundaries( trMesh );

for i = 1:numel(FF)
    % get the points of this boundary
    [coords] = freeBoundaryPoints(FF{i}, trMesh);
    
    % find the centre of the points
    centre = nanmean(coords, 1);
    
    % create a triRep of the boundary
    X = [centre; coords];
    numpts = size(X, 1);
    A = ones(numpts-1,1);
    B = [2:numpts]';
    C = [[3:numpts]' ; 2];
    TRI = [A B C];
    tr{i} = TriRep(TRI, X(:,1), X(:,2), X(:,3));
    
    % output values
    a(i) = sum(triarea(tr{i}));
    
    disp(['Perimeter is: ' num2str(l(i)) ' | Area is: ' num2str(a(i))]);
end

% plot
if plot
    drawFreeBoundaries(FF, trMesh);
    legend({'Boundary 1' 'Boundary 2' 'Boundary 3' 'Boundary 4' 'Boundary 5' 'Boundary 6' 'Boundary 7'})
    hold on
    hS = trisurf(getMesh(userdata));
    set(hS ...
        , 'facecolor', [.5 .5 .5] ...
        , 'edgecolor', 'none' ...
        , 'facealpha', .2 ...
        );
    set(gcf, 'color', 'w' ...
        );
    axis off equal vis3d;
end

end