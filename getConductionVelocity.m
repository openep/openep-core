function cvdata = getConductionVelocity( userdata )
% GETCONDUCTIONVELOCITY Returns the conduction velocity map of the 
% chamber
%
% Usage:
%   cvdata = getConductionVelocity( userdata )
% Where:
%   userdata  - see importcarto_mem
%   cvdata  - the conduction velocities, in m/s
%
% GETCONDUCTIONVELOCITY Calculate conduction velocities by calculating
% gradients of interpolated local activation times. GETCONDUCTIONVELOCITY
% makes use of a modified version of "Scattered Data Interpolation and 
% Approximation using Radial Base Functions" available from the Matlab
% FileExchange: Alex Chirokov (2020). Scattered Data Interpolation and 
% Approximation using Radial Base Functions 
% (https://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions), MATLAB Central File Exchange. Retrieved November 24, 2020.
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% cvdata = getConductionVelocity( userdata );
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

lats = userdata.electric.annotations.mapAnnot(getMappingPointsWithinWoI(userdata));
X = userdata.electric.egmX(getMappingPointsWithinWoI(userdata),:);
[~,~,~,~,cvdata] = RBFConductionVelocity(lats', X', getMesh(userdata).X');
cvdata = cvdata';

end