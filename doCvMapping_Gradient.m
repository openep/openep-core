function [cv, cvX, n, u] = doCvMapping_Gradient( userdata, int )
% DOCVMAPPING_GRADIENT Calculates conduction velocities using
% the gradient of the scalar activation time field
%
% Usage:
%   [cv, cvX, interpCv] = doCvMapping_Gradient( userdata, int )
% Where:
%   userdata - an OpenEP data structure
%   cv       - the calculated conduction velocity data, in m/s
%   cvX      - the Cartesian co-ordinates at which conduction velocity data
%              has been calculated. size(cvX) = [length(cv), 3].
%   u        - wave velocity vectors
%   n        - wave direction vectors (the unit vector field)
%
% DOCVMAPPING_GRADIENT Calculatess conduction velocities using the gradient
% of the interpolated local activation time field
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
d_interpLATs = grad';

% % calculate conduction velocities
% sgrad = grad .* grad;
% dp = sum(sgrad,2);
% cv = 1./sqrt(dp);

% calculate conduction velocities
mag_df=sqrt(d_interpLATs(1,:).^2+d_interpLATs(2,:).^2+d_interpLATs(3,:).^2);
cv=1./mag_df;
n=[cv.*d_interpLATs(1,:);cv.*d_interpLATs(2,:);cv.*d_interpLATs(3,:)];
u=[cv.*n(1,:);cv.*n(2,:);cv.*n(3,:)];

% transpose output
cv = cv';
n = n';
u = u';

end