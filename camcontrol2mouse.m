function camcontrol2mouse(h, varargin)
% CAMCONTROL2MOUSE enables mouse to control camera.
% Usage:
%   camcontrol2mouse(h)
%   camcontrol2mouse(h, 'off')
% Where:
%   h is a handle to axes or figure (all axes whithin will be used)
%
% CAMCONTROL2MOUSE let the user adjust the camera using the mouse. BUT all
% lights are held fixed relative to the camera. The ButtonDownFcn handles
% and related function handles of the axes are used.

% Author: Nick Linton (2011)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------
    
    if numel(h)>1
        for i = 1:numel(h)
            camcontrol2mouse(h(i), {varargin});
        end
        return
    end
                        
    type = get(h, 'Type');
    onOff = 'on';
    if nargin >= 2  &&  strcmpi(varargin{1},'off')
        onOff = 'off';
    end
    
    switch type
        case 'axes'
            set(h, 'ButtonDownFcn', @local_buttondown);
            hC = get(h, 'Children');
            for i = 1:numel(hC)
                camcontrol2mouse(hC(i), onOff);
            end
        case 'figure'
            hAx = get(h, 'Children');
            for i = 1:numel(hAx)
                camcontrol2mouse(hAx(i), onOff)
            end
        otherwise
            switch onOff
                case 'on'
                    set(h, 'ButtonDownFcn', @local_buttondown);
                case 'off'
                    set(h, 'ButtonDownFcn', '');
            end
    end
            
end

function local_buttondown(src, event)
    local_button(src, event, 'down')
end

function local_buttonup(src, event)
    local_button(src, event, 'up')
end

function local_button(src, event, keyword)
    persistent old_units
        
    %get the parent figure
    hFig = src;
    while ~isempty(hFig) && ~strcmpi('figure', get(hFig,'Type'))
        if strcmpi('axes', get(hFig,'Type'))
            hAx = hFig;
        end
        hFig = get(hFig,'Parent');
    end
    
    switch keyword
        case 'down'
            old_units = get(hFig, 'Units');
            set(hFig, 'Units', 'pixels');
            set(hFig, 'WindowButtonUpFcn', @local_buttonup)
            set(hFig, 'WindowButtonMotionFcn', @local_buttonmotion)
            local_buttonmotion(hFig, event, hAx)
        case 'up'
            set(hFig, 'WindowButtonMotionFcn','' , 'Units',old_units)
        otherwise
            error('local_button')
    end
end

function local_buttonmotion(src, event, varargin) 
    persistent hAx firstMousePoint hLights firstLightPos
    % debug 
    % persistent hTest
    
    hFig = src;
    
    cPosition = [0 0 0];	cTarget = [0 0 0];  cUpVector = [0 0 0];
    cViewVector = [0 0 0];  cRightVector = [0 0 0]; cUpVector = [0 0 0];
    cUpVector = [0 0 0];    cViewVector = [0 0 0];  cRightVector = [0 0 0];
        
    if nargin >= 3 %then this is the first call after button down
        hAx = varargin{1};
        firstMousePoint = get(hFig, 'CurrentPoint');
        hLights = findall(hAx, 'Type','light');
        firstLightPos = zeros(numel(hLights),3);
        for i = 1:numel(hLights)
            nested_getcamerageometry()
            
            pos = get(hLights(i), 'Position');
            c2Light = pos - cPosition;
            
            % express light position as RIGHT, UP, AWAY from camera
            firstLightPos(i,:) = [ sum(c2Light.*cRightVector) ...
                                    sum(c2Light.*cUpVector) ...
                                    sum(c2Light.*cViewVector) ...
                                   ];
            % debug
            % hTest = plot3dline([pos(:)';cPosition(:)'], 'r', 0.1, 5);
        end
    elseif isempty(firstMousePoint)
        return
    end
    point = get(hFig, 'CurrentPoint');
    
    delta = point-firstMousePoint;
    delta = -delta;
    
    dar = get(hAx, 'DataAspectRatio');
    nested_getcamerageometry()
    
    sel = get(src, 'SelectionType');
    if ~isempty(sel);   sel = lower(sel(1:3));  end
    switch sel
        case 'nor'
            [pos, up] = camrotate(cPosition,cTarget,dar,cUpVector, delta(1),delta(2), 'camera', 'y');
            set(hAx, 'cameraupvector', up, 'cameraposition', pos );
            % keep all lights in the same position relative to the camera
            nested_getcamerageometry();
            nested_adjustlights();
        case 'ext'
            delta = delta/120;
            delta(1) = -delta(1);
            [targ up] = camrotate(cTarget,cPosition,dar,cUpVector, delta(1),delta(2), 'camera', 'y');
            set(hAx, 'cameraupvector', up, 'cameratarget', targ );
            %              [lightPos, ~] = camrotate(pos,targ,dar,up,handles.userdata.view.camlightAz,handles.userdata.view.camlightEl,'camera','z');
            %              set(handles.userdata.hLight, 'position', lightPos, 'BusyAction', 'queue', 'Interruptible', 'on')
        case 'alt'
            k = 300;
            camzoom(hAx, (k-delta(2))/k );
    end
    firstMousePoint = point;
    
    function nested_getcamerageometry()
        cPosition = get(hAx, 'CameraPosition');
        cTarget = get(hAx, 'CameraTarget');
        cUpVector = get(hAx, 'CameraUpVector');
        
        cViewVector = cTarget - cPosition;
        cRightVector = cross(cViewVector, cUpVector);
        cUpVector = cross(cRightVector, cViewVector); % because the up vector isn't actually the up vector!!!
        
        % keep all vectors normalized
        cUpVector = cUpVector / norm2(cUpVector);
        cViewVector = cViewVector / norm2(cViewVector);
        cRightVector = cRightVector / norm2(cRightVector);
    end

    function nested_adjustlights()
            for iNested = 1:numel(hLights)
                nestedPos = firstLightPos(iNested,:);
                nestedPos = nestedPos(1)*cRightVector + nestedPos(2)*cUpVector + nestedPos(3)*cViewVector;
                nestedPos = nestedPos + cPosition;
                set(hLights(iNested), 'position', nestedPos, 'BusyAction', 'queue', 'Interruptible', 'on')
                % debug
                % delete(hTest); hTest = plot3dline([nestedPos(:)';cTarget(:)'], [1 0 0;0 1 0], 0.1, 5);
            end
    end
end
            
            
            
            