function [X, surfX] = getElectrogramX( userdata, varargin )
% GETELECTROGRAMX Returns the the electrode recording positions
%
% Usage:
%   C = getCentreOfMass( userdata, varargin )
% Where:
%   userdata   - see importcarto_mem
%   X - the 3D Cartesian co-ordinates
%   surfX - the surface-projected 3D Cartesian co-ordinates
%
% GETELECTROGRAMX accepts the following parameter-value pairs
%   'type'     {'bip'}|'uni'
%
% GETELECTROGRAMX can be used to access the Cartesian co-ordinates of the
% electrodes used to record the electrograms.
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% [X, surfX] = getElectrogramX(userdata, 'type', 'bip'); % default
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse input
nStandardArgs = 1; % UPDATE VALUE
type = 'bip';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'type'
                type = varargin{i+1};
        end
    end
end

switch lower(type)
    case 'bip'
        X = userdata.electric.egmX;
        surfX = userdata.electric.egmSurfX;
    case 'uni'
        if ~isfield(userdata.electric.egmUniX)
            error('OPENEP/GETELECTROGRAMX: There is no unipolar data associated with this data structure')
        else
            X = userdata.electric.egmUniX;
            surfX = findclosestvertex(getMesh(userdata), X);
        end
end