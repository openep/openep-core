classdef PointPicker < matlab.mixin.SetGet
    properties
        Mesh % the mesh A
        Surface % the plotted surface A
        hAx % the axes
        PointIndices % n x 2 array of fiducial indices, indexing into the mesh points lits
    end
    
    methods
        function obj = PointPicker(mesh, surface, ax)
            % Create a FiducialPicker for identifying pairs of fiducial markers on
            % pairs of surface meshes
            obj.Mesh = mesh;
            obj.hAx = ax;
            obj.Surface = surface;
            local_configurePicker(obj);
        end
        
        function keyPress(obj, src, event)
            disp('in key press')
        end
        
        function addPoint(obj, src, event)
                        
            set(src, 'Value', 0);
            
            DT(1) = get(obj.Surface, 'children');
            xyz = [DT(1).X DT(1).Y DT(1).Z];
            indexA = findclosestvertex(obj.Mesh, xyz);
            
            if isempty(obj.PointIndices)
                obj.PointIndices = [indexA];
            else
                obj.PointIndices(end+1, 1) = indexA;
            end
            
            axes(obj.hAx);
            obj.plotFiducialMarker(obj.Mesh, indexA)
            delete(DT)
        end
        
        function plotFiducialMarker(obj, m, ind)
            markerSize = 3;
            markerColor = 'r';
            hS = drawSphere(m.Points(ind,1), m.Points(ind,2), m.Points(ind,3), markerSize);
            set(hS, 'facecolor', markerColor);           
        end
        
        function local_configurePicker(obj)
            tb = axtoolbar(gca, {'rotate','datacursor'});
            btn = axtoolbarbtn(tb, 'state');
            btn.Icon = 'mygridicon.png';
            btn.Tooltip = 'Add Point';
            btn.ValueChangedFcn = @obj.addPoint;
            hold on
        end
        
        function saveState(obj)
            fidInds = obj.PointIndices;
            save('fiducialpickerstate' ...
                , 'fidInds');                
        end
        
        function loadState(obj)
            load('pointpickerstate');
            obj.PointIndices = fidInds;
            axes(obj.hAx);
            for i = 1:size(obj.PointIndices,1)
                obj.plotFiducialMarker(obj.Mesh, obj.PointIndices(i,1));
            end
        end
    end
end
