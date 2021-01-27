function [tr2, isVertUsed] = editTriangulation(tr)
% EDITTRIANGULATION Graphically remove triangles from a TriRep object
%
% Usage:
%   tr2 = editTriangulation(tr)
% Where:
%   tr  - is the original triangulation
%   tr2 - is the new triangulation with elements removed
%   isVertUsed - indexes into tr.X for vertices that are used in the new
%   triangulation, tr2
%
% EDITTRIANGULATION uses GET3DFACES to remove triangles from a TriRep
% object. Controls:
%   Left click          - select triangles to remove
%   Shift-Left click    - select triangles to keep
%   Ctrl-Left click     - select area up to the boundary
%   d                   - done
%
% Author: Steven Williams (2013)
% Modifications -
%   2016 - refactored by moving userdata manipulation into a separate fcn
%
% Info on Code Testing:
% ---------------------
%
% ---------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% Edit to remove the valve orifice
disp('E D I T   T R I A N G U L A T I O N');
disp('-----------------------------------');
disp('Left click - select triangles to remove (shift to undo)');
disp('Cltr-left - select area');
disp('Press d when done');
disp('');

editFig = figure;

% Construct the trirep
no = tr.X; %node
el = tr.Triangulation; %element

hP = trisurf(tr);
set(hP, 'FaceVertexCData', zeros(length(tr.Triangulation),3)+0.7 ...
    ,'FaceColor', 'flat' ...
    );
axis equal vis3d
light

get3dfaces(hP,'on');

%Wait until a key is pressed, then only progressif it was 'd'
key = get(editFig, 'CurrentCharacter');
while ~strcmpi(key, 'd')
    pause;
    key = get(editFig, 'CurrentCharacter');
end

% Now, continue by creating the new TriRep
iRemove = get3dfaces(hP);
close(editFig);

el(iRemove, :) = [];
tr2 = repack(TriRep(el, no));

[tr2, isVertUsed] = repack(tr2);
end