function [cv, cvX, u, n] = getConductionVelocity( userdata, varargin )
% GETCONDUCTIONVELOCITY Returns the conduction velocity map of the chamber
%
% Usage:
%   cvdata = getConductionVelocity( userdata )
%   cvdata = getConductionVelocity( userrdata, interpolator)
%   cvdata = getConductionVelocity( ... , param, value )
% Where:
%   userdata - an OpenEP data structure
%   cv       - the calculated conduction velocity data, in m/s
%   cvX      - the Cartesian co-ordinates at which conduction velocity data
%              has been calculated. size(cvX) = [length(cv), 3].
%   u        - wave velocity vectors
%   n        - wave direction vectors (the unit vector field)
%
% GETCONDUCTIONVELOCITY accepts the following parameter-value pairs
%   `'method'`        `'triangulation'` | `'cosinefit'` | `{'radialbasis'}` |
%                     `'gradient'`
%   `'interpolator'`  `{'scatteredInterpolant'}` | `'radialbasis'` |
%                     `'localsmoothing'`; or an instance of openEpDataInterpolator
%   TODO add the rbf options
%
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

CVLIMIT = 10; %m/s
DISTANCETHRESHOLD = 10; %mm

nStandardArgs = 1; % UPDATE VALUE
method = 'radialbasis';
interpolator = 'scatteredinterpolant'; %COS: Changed default from scatteredinterpolant
rbfoptions.doOptimisation=false; %defaults
rbfoptions.shapeParameter=1;
rbfoptions.basisFunction='multiquadratic';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'method'
                method = varargin{i+1};
            case 'interpolator'
                interpolator = varargin{i+1};
            case 'dooptimisation'
                doOptimisation = varargin{i+1};
                rbfoptions.doOptimsation=doOptimisation;
            case 'shapeparameter'
                shapeParameter = varargin{i+1};
                rbfoptions.shapeParameter=shapeParameter;
            case 'basisfunction'
                basisFunction = varargin{i+1};
                rbfoptions.basisFunction=basisFunction;
            case 'distancethreshold'
                DISTANCETHRESHOLD = varargin{i+1};
            case 'cvlimit'
                CVLIMIT = varargin{i+1};
        
        end
    end
end

% first create an interpolator if needed; only for gradient method
if strcmpi(method, 'gradient')
    if ischar(interpolator)
        int = openEpDataInterpolator(interpolator);
    else
        int = interpolator;
    end
end
int.distanceThreshold = DISTANCETHRESHOLD;

% %RBF option setting  - COS
% rbfoptions.doOptimisation=false; %defaults
% rbfoptions.shapeParameter=1;
% rbfoptions.basisFunction='multiquadratic';
% if nargin > nStandardArgs
%     if ischar(varargin{2})
%         for i = 1:2:nargin-nStandardArgs
%             switch varargin{i}
%                 case 'doOptimisation'
%                     doOptimisation = varargin{i+1};
%                     rbfoptions.doOptimsation=doOptimisation;
%                 case 'shapeParameter'
%                     shapeParameter = varargin{i+1};
%                     rbfoptions.shapeParameter=shapeParameter;
%                 case 'basisFunction'
%                     basisFunction = varargin{i+1};
%                     rbfoptions.basisFunction=basisFunction;
%             end
%         end
%     end
% end

switch lower(method)
    case 'triangulation'
        [cv, cvX, u, n] = doCvMapping_Triangulation(userdata);

    case 'cosinefit'
        [cv, cvX, u, n] = doCvMapping_CosineFit(userdata);

    case 'radialbasis'
        [cv, cvX, u, n] = doCvMapping_RadialBasis(userdata, rbfoptions);

    case 'gradient'
        [cv, cvX, u, n] = doCvMapping_Gradient(userdata, int);
end

% accept only those conduction velocity values in proximity to electrodes
vtx = getVerticesNearMappingPoints(userdata, DISTANCETHRESHOLD);
cv(~vtx) = NaN;
cvX(~vtx,:) = NaN;
n(~vtx,:) = NaN;
u(~vtx,:) = NaN;
disp(['OPENEP/GETCONDUCTIONVELOCITY: ' num2str(sum(~vtx)) ' CV values were removed which were more than ' num2str(DISTANCETHRESHOLD) 'mm from a mapping point']);

% remove any non physiological values over the CVLIMIT
isOverCvLimit = cv>CVLIMIT;
cv(isOverCvLimit) = NaN;
cvX(isOverCvLimit,:) = NaN;
n(isOverCvLimit,:) = NaN;
u(isOverCvLimit,:) = NaN;
disp(['OPENEP/GETCONDUCTIONVELOCITY: ' num2str(sum(isOverCvLimit)) ' CV values were removed which were greater than ' num2str(CVLIMIT) 'm/s']);


end