function [cv, cvX, interpCv] = doCvMapping_Gradient( userdata, int )
% DOCVMAPPING_GRADIENT Calculates conduction velocities using
% the gradient of the scalar activation time field
%
% Usage:
%   [cv, cvX, interpCv] = doCvMapping_Gradient( userdata, int )
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
% DOCVMAPPING_GRADIENT Calculatess conduction velocities using
% omnipoles
%
% Author: 
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% -----------------------------------s----------------------------
% load openep_dataset_1.mat
% int = openEpDataInterpolator
% [cv, cvX, interpCv] = getConductionVelocity(userdata, 'method', 'gradient', 'interpolator', int);
% histogram(cv)
% drawMap(userdata, 'data', interpCv, 'type', 'cv')
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

% perform global interpolation of the activation field
latMap = int.interpolate(X, lats, getVertices(userdata));

% calculate the gradeint of this field
[cvX, grad] = trigrad(getMesh(userdata), latMap);

% calculate conduction velocities
sgrad = grad .* grad;
dp = sum(sgrad,2);
cv = 1./sqrt(dp);

% accept only those conduction velocity values in proximity to electrodes
vtx = getVerticesNearMappingPoints(userdata, DISTANCETHRESHOLD);
cv(~vtx) = [];
cvX(~vtx,:) = [];
disp(['OPENEP/DOCVMAPPING_GRADIENT: ' num2str(sum(~vtx)) ' CV values were removed which were more than ' num2str(DISTANCETHRESHOLD) 'mm from a mapping point']);

% remove any non physiological values over the CVLIMIT
isOverCvLimit = cv>CVLIMIT;
cv(isOverCvLimit) = [];
cvX(isOverCvLimit,:) = [];
disp(['OPENEP/DOCVMAPPING_GRADIENT: ' num2str(sum(isOverCvLimit)) ' CV values were removed which were greater than ' num2str(CVLIMIT) 'm/s']);

% interpolate
interpCv = int.interpolate(cvX, cv, getVertices(userdata));

end