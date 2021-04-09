function int = createInterpolator( type, options )
% CREATEINTERPOLATOR Creates an interpolation scheme for use with OpenEP
% data
%
% Usage:
%   int = createInterpolator( type )
%
% Where:
%   type - the type of interpolator required which may be one of, 
%          `'scatteredInterpolant'`, `'radialbasis'` or `'localsmoothing'`
%   int - a function handle to the interpolator
%   options - a structure with the following fields:
%       options.method
%       options.extrapolationmethod
%       options.distancethreshold
%
% CREATEINTERPOLATOR creates and interpolation scheme for use with OpenEP
% data
%
% Author: Steven Williams (2021) (Copyright)
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

switch(type)
    
    case 'scatteredInterpolant'
        
        method = options.method;
        extrapolationMethod = options.extrapolationMethod;
        distanceThreshold = options.distanceThreshold;
        
    case 'radialbasis'
        
    case 'localsmoothing'
        
end