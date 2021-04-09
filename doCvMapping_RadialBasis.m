function [cv, cvX, interpCv] = doCvMapping_RadialBasis( userdata, int )
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

% first perform global interpolation using radial basis function to 
% calculate conduction velocities
lats = userdata.electric.annotations.mapAnnot(getMappingPointsWithinWoI(userdata));
X = userdata.electric.egmX(getMappingPointsWithinWoI(userdata),:);
cvX = userdata.surface.triRep.X;
[~,~,~,~,cvdata] = RBFConductionVelocity(lats', X', cvX');
cv = cvdata';

% accept only those conduction velocity values in proximity to electrodes?
% TODO. We should in some way limit cv and cvX only to values that are
% likely to be real; i.e. in close proximity; or at; mapping points.

% now do interpolation, using the interpolator specified, so we have a full
% dataset at each of the mesh nodes.
vtx = getVertices(userdata, 'used', false);
interpCv = int.interpolate(X, cv, vtx);

end