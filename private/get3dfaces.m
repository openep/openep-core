function iFaces = get3dfaces(varargin)
% GET3DFACES lets the user select faces from a patch object interactively.
% Usage:
%   iFaces = get3dfaces(hPatch)
%   iFaces = get3dfaces(hPatch, 'on')
%   iFaces = get3dfaces(hPatch, 'off')
%   iFaces = get3dfaces(hPatch, 'override', selectedfaces)
% Where:
%   iFaces is a list of the faces that have already been selected
%   'on' and 'off' turn the utility on/off

% Author: Nick Linton (2011)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------


% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------
    cSelect = [0 0 1];
    
    hPatch = varargin{1};
    hFig = get(hPatch, 'Parent');
    while ~strcmpi(get(hFig, 'Type'), 'figure')
        hFig = get(hFig, 'Parent');
    end
    
    tag = get(hPatch, 'Tag');
    if nargin >= 2
        switch varargin{2}
            case 'on'
                turn = 'on';
            case 'off'
                turn = 'off';
            case 'override'
                userdata = get(hPatch, 'Userdata');
                userdata.iFacesSelected = varargin{3};
                
                faceColor = userdata.originalFaceColor;
                faceColor(userdata.iFacesSelected, :) = repmat(cSelect, length(userdata.iFacesSelected), 1);
                
                set(hPatch, 'FaceVertexCData', faceColor)
                set(hPatch, 'Userdata', userdata);
                
                turn = 'on';
        end
    elseif strcmp(tag,'get3dface')
        turn = 'off';
    else
        turn = 'on';
    end
    
    %change the colors to direct TrueColor
    faceColor = get(hPatch, 'FaceVertexCData');
    faces = get(hPatch, 'Faces');
    if length(faceColor)~=length(faces) || min(size(faceColor)) == 1
        warning('GET3DFACES: only works on patches shaded by face with TrueColor.')
        set(hPatch ...
                ,'FaceVertexCData', zeros(length(faces),3)+0.5 ...
                ,'FaceColor', 'flat' ...
                                        )
    end
        

    if strcmp(turn,'off')
        set(hFig, 'WindowButtonDownFcn', '')
        set(hFig,'Pointer','arrow')
        set(hPatch, 'Tag', '')
    else
        set(hFig, 'WindowButtonDownFcn', @buttondown)
        set(hFig,'Pointer','circle')
        set(hPatch, 'Tag', 'get3dface')
    end
    userdata = get(hPatch,'Userdata');
    if ~isempty(userdata)
        iFaces = userdata.iFacesSelected;
    else
        iFaces = [];
    end
end



function buttondown(src,evnt) %#ok<*INUSD>
    set(src, 'WindowButtonUpFcn', @buttonup)
    set(src, 'WindowButtonMotionFcn', @selectpoint)
    set(src,'Pointer','circle')
    selectpoint(src,[])
end

function buttonup(src,evnt)
    set(src, 'WindowButtonMotionFcn', '')
    set(src, 'WindowButtonUpFcn', '')
end


function selectpoint(src,evnt)
    if ~strcmpi(get(gco, 'Tag'), 'get3dface');  return; end
    hold on
    
    cSelect = [0 0 1];
    
    userdata = get(gco,'Userdata');
    if isempty(userdata)
        userdata.iFacesSelected = [];
        userdata.originalFaceColor = get(gco, 'FaceVertexCData');
    end

    [~, ~, ~, ~, faceIndex] = select3d();
    
    if strcmpi(get(src,'SelectionType'), 'normal')
        userdata.iFacesSelected = unique([userdata.iFacesSelected ; faceIndex]);
    elseif strcmpi(get(src,'SelectionType'), 'extend')
        userdata.iFacesSelected(userdata.iFacesSelected == faceIndex) = [];
    elseif strcmpi(get(src,'SelectionType'), 'alt')
        % now we need to spread to include all the triangles, without
        % crossing previously selected triangles
        if any(faceIndex==userdata.iFacesSelected)
            return
        end
        faces = get(gco, 'Faces');
        nFaces = length(faces);
    

        iFacesSelected = userdata.iFacesSelected;
        iFacesNotSelected = 1:nFaces;
        iFacesNotSelected(iFacesSelected) = [];
        
        facesNotSelected = faces(iFacesNotSelected,:);
        start = find(faceIndex == iFacesNotSelected, 1, 'first');
        [~, newFaceIndices]  = getallattachedsurface(facesNotSelected,  start, 'face');
        
        userdata.iFacesSelected = unique([userdata.iFacesSelected ; iFacesNotSelected(newFaceIndices)']);
    end
    
    faceColor = userdata.originalFaceColor;
    faceColor(userdata.iFacesSelected, :) = repmat(cSelect, length(userdata.iFacesSelected), 1);
    
    set(gco, 'FaceVertexCData', faceColor)
    set(gco,'Userdata', userdata);
    
end