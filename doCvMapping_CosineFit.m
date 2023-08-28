function [cv, cvX, interpCv] = doCvMapping_CosineFit( userdata, int )
% DOCVMAPPING_COSINEFIT Calculates conduction velocities using
% triangulation of electrode activation times
%
% Usage:
%   [cv, cvX, interpCv] = doCvMapping_CosineFit( userdata, int )
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
% DOCVMAPPING_COSINEFIT Calculatess conduction velocities using
% cosine fit technique
%
% Author: 
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% -----------------------------------s----------------------------
%
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% first perform global interpolation using cosine fit to 
% calculate conduction velocities

% TODO

% accept only those conduction velocity values in proximity to electrodes?
% TODO. We should in some way limit cv and cvX only to values that are
% likely to be real; i.e. in close proximity; or at; mapping points.

% now do interpolation, using the interpolator specified, so we have a full
% dataset at each of the mesh nodes.
vtx = getVertices(userdata, 'used', false);
interpCv = int.interpolate(X, cv, vtx);

end