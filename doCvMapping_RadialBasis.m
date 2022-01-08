function [cv, cvX, u, n] = doCvMapping_RadialBasis( userdata, rbfoptions)
% DOCVMAPPING_RADIALBASIS Calculates conduction velocities using radial
% basis functions
% Usage:
%   [cvX, cv, u, n] = doCvMapping_RadialBasis( userdata, rbfoptions )
% Where:
%   userdata - an OpenEP data structure
%   cv       - the calculated conduction velocity data, in m/s
%   cvX      - the Cartesian co-ordinates at which conduction velocity data
%              has been calculated. size(cvX) = [length(cv), 3].
%   u        - wave velocity vectors
%   n        - wave direction vectors (the unit vector field)
%   rbfoptions  - TODO COMPLETE THE DOCUMENTATION FOR this structure
%   
% DOCVMAPPING_RADIALBASIS Calculates conduction velocities using radial
% basis functions
%
% Author: Steven Willians / Chris O'Shea (2022) (Copyright)
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

% get relevant mapping points
X = userdata.electric.egmX(getMappingPointsWithinWoI(userdata),:);

% calculate the local activation times at each mapping point
mapAnnot = userdata.electric.annotations.mapAnnot(getMappingPointsWithinWoI(userdata));
refAnnot = userdata.electric.annotations.referenceAnnot(getMappingPointsWithinWoI(userdata));
lats = mapAnnot - refAnnot;

%get all mesh points
x_all=getVertices(userdata);

% perform global interpolation of the activation field
% latMap = int.interpolate(X, lats, x_all);
% size(latMap)

% fit rbf field to the lat map
[interpLATs, d_interpLATs]=rbf_interpolator(lats,X,x_all,rbfoptions);

% calculate the gradeint of this field
%[cvX, d_interp_values] = trigrad(getMesh(userdata), interp_values);

% calculate conduction velocities
cvX = x_all;
mag_df=sqrt(d_interpLATs(1,:).^2+d_interpLATs(2,:).^2+d_interpLATs(3,:).^2);
cv=1./mag_df;
n=[cv.*d_interpLATs(1,:);cv.*d_interpLATs(2,:);cv.*d_interpLATs(3,:)];
u=[cv.*n(1,:);cv.*n(2,:);cv.*n(3,:)];

% transpose output
cv = cv';
n = n';
u = u';

% calculate conduction velocities
% sgrad = d_interp_values .* d_interp_values;
% dp = sum(sgrad,2);
% cv = 1./sqrt(dp);

% We decided no longer to return the interpolated field
% % now do interpolation, using the interpolator specified, so we have a full
% % dataset at each of the mesh nodes.
% interpCv = int.interpolate(cvX, cv, getVertices(userdata));

end