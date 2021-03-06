function tr = getClosedSurface( userdata )
% GETCLOSEDSURFACE Fills all the holes in the userdata surface
%
% Usage:
%   tr = getClosedSurface( userdata )
% Where:
%   userdata  - see importcarto_mem
%   tr  - a triRep object
%
% GETCLOSEDSURFACE Returns a new surface representation of the anatomical 
% model with all the holes in the mesh filed. Closes the surface by the 
% following algorithm. First, every complete free boundary is identified. 
% Second, the barycentre of the free boundary is identified. Third, a 
% triangulation is created covering this hole. Finally, the additional 
% triangles are added to the TriRep.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% tr = getClosedSurface( userdata )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% get the trirep
tr = userdata.surface.triRep;

% get all the free boundaries
[ FF, ~ ] = freeBoundaries( tr );

% step through each free boundary
for i = 1:numel(FF)
    % get the points of this boundary
    coords = freeBoundaryPoints(FF{i}, tr);
    
    % find the centre of these points
    centre = nanmean(coords, 1);
    
    % create a triRep of the hole
    X = [centre; coords];
    numpts = size(X, 1);
    A = ones(numpts-1,1);
    B = [2:numpts]';
    C = [[3:numpts]' ; 2];
    TRI = [A B C];
    newtr = TriRep(TRI, X(:,1), X(:,2), X(:,3));
    
    % add this new triRep to the origianl triRep
    tr = addtrirep(tr, newtr);
end

end