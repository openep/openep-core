function [cv, cvX, interpCv] = getConductionVelocity( userdata, varargin )
% GETCONDUCTIONVELOCITY Returns the conduction velocity map of the chamber
%
% Usage:
%   cvdata = getConductionVelocity( userdata )
%   cvdata = getConductionVelocity( userrdata, interpolator)
%   cvdata = getConductionVelocity( ... , param, value )
% Where:
%   userdata - see importcarto_mem
%   cv       - the calculated conduction velocity data, in m/s
%   cvX      - the Cartesian co-ordinates at which conduction velocity data
%              has been calculated. size(cvX) = [length(cv), 3].
%   interpCv - conduction velocity data interpolated across the surface of
%              the shell.
%              size(interpCv) = [length(userdata.surface.triRep.X), 1].
%
% GETCONDUCTIONVELOCITY accepts the following parameter-value pairs
%   `'method'`        `'triangulation'` | `'cosinefit'` | `{'radialbasis'}` |
%                     `'omnipole'` | `'gradient'` | `'eikonal'`
%   `'interpolator'`  `{'scatteredInterpolant'}` | `'radialbasis'` |
%                     `'localsmoothing'`

% GETCONDUCTIONVELOCITY Returns the conduction velocity data of the chamber.
% Five methods for calculating conduction velocity are provided as
% described below. The conduction velocity data interpolated over the
% surface of the shell is also returned. The desired interpolation scheme
% can be set by passing in an OpenEPInterpolator as the second argument.
%
% Conduction velocity methods
% ---------------------------
% (*) triangulation
%   description goes here
%
% (*) cosinefit
%   description goes here
%
% (*) radialbasis
%   description goes here
%
% (*) omnipole
%   description goes here
%
% (*) gradient
%   description goes here
%
% (*) eikonal
%   description goes here
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

nStandardArgs = 1; % UPDATE VALUE
method = 'radialbasis';
interpolator = 'scatteredinterpolant';
if nargin > nStandardArgs
    if ischar(varargin{2})
        for i = 1:2:nargin-nStandardArgs
            switch varargin{i}
                case 'method'
                    method = varargin{i+1};
                case 'interpolator'
                    interpolator = varargin{i+1};
            end
        end
    else
        interpolator = varargin{2};
    end
end

% first create an interpolator
if ischar(interpolator)
    int = createInterpolator(interpolator);
else
    int = interpolator;
end
    
    switch lower(method)
        case 'triangulation'
            [cv, cvX, interpCv] = doCvMapping_Triangulation(userdata, int);
            
        case 'cosinefit'
            [cv, cvX, interpCv] = doCvMapping_CosineFit(userdata);
            
        case 'radialbasis'
            [cv, cvX, interpCv] = doCvMapping_RadialBasis(userdata);
            
            %         lats = userdata.electric.annotations.mapAnnot(getMappingPointsWithinWoI(userdata));
            %         X = userdata.electric.egmX(getMappingPointsWithinWoI(userdata),:);
            %         [~,~,~,~,cvdata] = RBFConductionVelocity(lats', X', userdata.surface.triRep.X');
            %         cvdata = cvdata';
            
        case 'omnipole'
            [cv, cvX, interpCv] = doCvMapping_Omnipole(userdata, int);
            
        case 'gradient'
            [cv, cvX, interpCv] = doCvMapping_Gradient(userdata, int);
            
        case 'eikonal'
            [cv, cvX, interpCv] = doCvMapping_Eikonal(userdata, int);
            
    end
    
end