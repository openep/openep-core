function hL = cameraLight
% CAMERALIGHT Adds light to the current axis.
% Usage:
%   cameralight
%   cameralight(param1, Value1, ..., ParamN, ValueN)
%   L = cameralight(...)
%
% CAMERALIGHT Positions a light object at the camera position and listens
% to update the position whenever the camera moves. For parameter-value 
% pairs see light.
%
% Author: Steven Williams (2014)
% Modifications -
%   Changed listener to 'MarkedClean' after R2014b removed the ability to
%   listen post-set to handle graphic object properties
%       - Steven Williams (2016)
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

hAx = gca();

hL = camlight();

positionLight();

el_global = addlistener(hAx, 'MarkedClean', @positionLight);

    function positionLight(varargin)
        campos = get(hAx, 'cameraposition');
        try
            set(hL, 'position', campos);
        catch
            disp('About to delete the listener');
            delete(el_global);
        end
    end
end