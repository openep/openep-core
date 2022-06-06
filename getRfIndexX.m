function X = getRfIndexX( userdata, iLesion, varargin )
% GETRFINDEXX Returns the Cartesian co-ordinates of an RF ablation tag
%
% Usage:
%   X = getRfIndexX( userdata, iLesion, varargin )
% Where:
%   userdata  - an OpenEP data structure
%   iLesion  - the index of the ablation lesion
%   X - Cartesion co-ordinates of the lesion
%
% GETRFINDEXX accepts the following parameter-value pairs
%   'type'      {'3d'} | 'surface'
%
% GETRFINDEXX Returns the Cartesian co-ordinates of the ablation lesion
% specified in usredata.rfindex, with the index of iLesion. By default the
% 3d (actual) locations of the lesion are returned. If required the
% co-ordinates projected to the surface can be returned instead.
%
% Author: Steven Williams (2022) (Copyright)
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

nStandardArgs = 2; % UPDATE VALUE
type = '3d';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'type'
                type = varargin{i+1};
        end
    end
end

switch lower(type)
    case '3d'
        X = userdata.rfindex.tag.X(iLesion,:);

    case 'surface'
        X = userdata.rfindex.tag.X(iLesion,:);

        % find the surface projection
        warning('off', 'FINDCLOSESTVERTEX:hasToTrim')
        [closestVertices, ~] = findclosestvertex(getMesh(userdata), X, true);
        warning('on', 'FINDCLOSESTVERTEX:hasToTrim')

        % Now work out the surface projections
        allVertices = getVertices(userdata);
        X = allVertices(closestVertices,:);
end

