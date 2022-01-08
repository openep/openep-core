function h = plotQuiverField(userdata)
% PLOTQUIVERFIELD Brief description goes here
%
% Usage:
%   h = plotQuiverField(userdata)
% Where:
%   userdata - an OpenEP data structure
%   h        -  a quivergroup handle
%
% PLOTQUIVERFIELD accepts the following parameter-value pairs
%   'rbfoptions'     see doCvMapping_RadialBasis
%
% PLOTQUIVERFIELD plots a quiver field of activation directions into the
% current axis
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

% calculate the conduction vectors
[cv, cvX, n, u] = getConductionVelocity(userdata, 'method', 'radialbasis');

% get the surface normals
[normals, userdata] = getNormals(userdata);
iV = findclosestvertex(getMesh(userdata, 'trirep'), cvX);
normals = normals(iV,:);

% project the velocities in tangent direction
projU=projectVertexVectors(cvX,u,normals);

% compute unit-vectors:
uv=createUnitVectors(projU);

% plot
hold on
h = quiver3(cvX(:,1),cvX(:,2),cvX(:,3),uv(:,1),uv(:,2),uv(:,3) ...
    , 'color', 'k' ...
    , 'linewidth', 1 ...
    );

end