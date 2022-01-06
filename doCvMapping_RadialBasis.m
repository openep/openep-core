function [cv, cvX, interpCv] = doCvMapping_RadialBasis( userdata, int, rbfoptions)
% DOCVMAPPING_RADIALBASIS Calculates conduction velocities using radial
% basis functions
%
% Usage:
%   [cv, cvX, interpCv] = doCvMapping_RadialBasis( userdata )
% Where:
%   userdata - see importcarto_mem
%   int      - see openEpDataInterpolator.m
%   cv       - the calculated conduction velocity data, in m/s
%   cvX      - the Cartesian co-ordinates at which conduction velocity data
%              has been calculated. size(cvX) = [length(cv), 3].
%   interpCv - conduction velocity data interpolated across the surface of
%              the shell.
%              size(interpCv) = [length(userdata.surface.triRep.X), 1].
%   
% DOCVMAPPING_RADIALBASIS Calculates conduction velocities using radial
% basis functions
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

CVLIMIT = 10; %m/s
DISTANCETHRESHOLD = 10; %mm

% get relevant mapping points
X = userdata.electric.egmX(getMappingPointsWithinWoI(userdata),:);

% calculate the local activation times at each mapping point
mapAnnot = userdata.electric.annotations.mapAnnot(getMappingPointsWithinWoI(userdata));
refAnnot = userdata.electric.annotations.referenceAnnot(getMappingPointsWithinWoI(userdata));
lats = mapAnnot - refAnnot;

%get all mesh points
x_all=getVertices(userdata);

% perform global interpolation of the activation field
latMap = int.interpolate(X, lats, x_all);
size(latMap)
% fit rbf field to the lat map
[interpLATs]=rbf_interpolator(latMap,x_all,x_all,rbfoptions);

% calculate the gradeint of this field
[cvX, grad] = trigrad(getMesh(userdata), interpLATs);

% calculate conduction velocities
sgrad = grad .* grad;
dp = sum(sgrad,2);
cv = 1./sqrt(dp);

% accept only those conduction velocity values in proximity to electrodes
vtx = getVerticesNearMappingPoints(userdata, DISTANCETHRESHOLD);
cv(~vtx) = [];
cvX(~vtx,:) = [];
disp(['OPENEP/DOCVMAPPING_RADIALBASIS: ' num2str(sum(~vtx)) ' CV values were removed which were more than ' num2str(DISTANCETHRESHOLD) 'mm from a mapping point']);

% remove any non physiological values over the CVLIMIT
isOverCvLimit = cv>CVLIMIT;
cv(isOverCvLimit) = [];
cvX(isOverCvLimit,:) = [];
disp(['OPENEP/DOCVMAPPING_RADIALBASIS: ' num2str(sum(isOverCvLimit)) ' CV values were removed which were greater than ' num2str(CVLIMIT) 'm/s']);

% now do interpolation, using the interpolator specified, so we have a full
% dataset at each of the mesh nodes.
interpCv = int.interpolate(cvX, cv, getVertices(userdata));

end