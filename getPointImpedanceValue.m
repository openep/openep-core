function a = getPointImpedanceValue(imp, tim)
% GETPOINTIMPEDANCEVALUE Provides an algorithm for giving point impedance
% values
%
% Usage:
%   h = myfunction(b)
% Where:
%   a - is the output
%   imp - impedance values at times tim
%
% Impedance values are streaming every 100ms to Carto® 3 system from the 
% RF Generator. A time range of -7.5s to +2.5s is output for each point and
% saved in userdata.electric.impedances (.time and .value) in the mat
% files.
%
% GETPOINTIMPEDANCEVALUE converts these time series to a single value by
% taking the peak impedance immediately prior to the 0 time point.
%
% Author: Steven Williams (2014) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
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


[~,loc] = findpeaks(imp);
A = [tim(loc) imp(loc)];
A(A(:,1)>0,:) = [];
a = A(end,2);

end

