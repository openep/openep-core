classdef HoleCutter < matlab.mixin.SetGet
    properties
        Mesh % the mesh
        Surface % the plotted surface
        DataTips % the datatips
        Method = 'vertex'; % whether to use 'vertex' or 'centroid' of the datatips
        XYZ % co-ordinates of datatips
        iV % index of the nodes corresponding to DataTips
        centroidXYZ % co-ordinates of the closest centroid to each datatip
        iCentroid % index of each centroid closest to each datatip
        Route % the route around the DataTips across the mesh
        eCross % edges crossing the Route
        Algorithm % the geodesic algorithm
        Path % cell array of paths between adjacent points created by exact_geodesic library
        PathX % Cartesian co-ordinates of the path
        PathNodeDistances % nx2 array of the closests nodes and the distances to those nodes
        PathHandles % array of handles to any plotted paths
    end
    
    methods
        function obj = HoleCutter(trMesh, hT)
% Create a HoleCutter object for cutting holes in a mesh
% Begin by using the DataTip function with <shift><leftclick> to select n>4
% points on the mesh. Then instantiate a HoleCutter object:
%             hc = HoleCutter(trMesh2, hT);
%
% Then call the following functions:
%             hc.getDataTips();
%             hc.getDataTipCoords();
%             hc.initialiseGeodesicCalculation();
%             hc.computePointOrder();
%             hc.computePath();
%             hc.getPathX();
%             hc.plotPath();
%             hc.deleteGeodesic();
% 
% If needed, data tips can be added or moved. When data tips are moved, the
% geodesic path will re-calculate.
            
            obj.Mesh = trMesh;
            obj.Surface = hT;
            set(obj.Surface, 'ButtonDownFcn', @obj.valueChangedFcn);
        end
        
        function getDataTips(obj)
            % Get the datatips from the surface, store them in the
            % HoleCutter object and set their ValueChangedFcn
            dataTips = get(obj.Surface, 'children');
            obj.DataTips = dataTips;
            set(obj.DataTips, 'ValueChangedFcn', @obj.valueChangedFcn);
        end
        
        function initialiseGeodesicCalculation(obj)
            mesh = geodesic_new_mesh(obj.Mesh.Points, obj.Mesh.ConnectivityList);
            obj.Algorithm = geodesic_new_algorithm(mesh, 'exact');
        end
        
        function computePointOrder(obj)
            % Order the points
            numPts = numel(obj.iV);
            route = NaN(numPts,2);
            route(1) = 1;
            route(end) = 1;
            for i = 1:numPts-1
                clear source_points;
                source_points{1} = geodesic_create_surface_point('vertex', obj.iV(route(i,1)), obj.Mesh.Points(obj.iV(route(i,1)),:));
                geodesic_propagate(obj.Algorithm, source_points);
                [~, distances] = geodesic_distance_and_source(obj.Algorithm);
                D = distances(obj.iV);
                toRemove = unique(route);
                toRemove(isnan(toRemove)) = [];
                D(toRemove) = Inf;
                [~,iNextPoint] = min(D);
                route(i,2) = iNextPoint;
                route(i+1,1) = iNextPoint;
            end
            obj.Route = route;
        end
        
        function labelDataTips(obj)
            % Label the DataTips (not currently used)
            for i = 1:numel(obj.DataTips)
                obj.DataTips(i).UserData = i;
            end
        end
        
        function getDataTipCoords(obj)
            % Get the co-ordinates of the data tips and identify the
            % closest mesh vertex
            for i = 1:numel(obj.DataTips)
                obj.XYZ(i,1:3) = [obj.DataTips(i).X obj.DataTips(i).Y obj.DataTips(i).Z];
                obj.iV(i) = findclosestvertex(obj.Mesh, obj.XYZ(i,:));
            end
        end
        
        function convertDataTipCoordsToCentroids(obj)
            % Convert datatips to centroid locations
            [~, allCentroids] = tricentroid(obj.Mesh);
            for i = 1:numel(obj.DataTips)
                obj.iCentroid(i) = findclosestvertex(allCentroids, obj.XYZ(i,1:3));
                obj.centroidXYZ(i,1:3) = allCentroids(obj.iCentroid(i),:);
            end
        end
         
        function computePath(obj)
            % Compute the geodesic path between all the datatips
            clear source_points
            switch obj.Method
                case 'vertex'
                    for i = 1:size(obj.Route,1)
                        source_points{1} = geodesic_create_surface_point('vertex', obj.iV(obj.Route(i,1)), obj.Mesh.Points(obj.iV(obj.Route(i,1)),:));
                        geodesic_propagate(obj.Algorithm, source_points);
                        destination = geodesic_create_surface_point('vertex', obj.iV(obj.Route(i,2)), obj.Mesh.Points(obj.iV(obj.Route(i,2)),:));
                        obj.Path{i} = geodesic_trace_back(obj.Algorithm, destination);
                    end
                case 'centroid'
                    for i = 1:size(obj.Route,1)
                        source_points{1} = geodesic_create_surface_point('face',obj.iCentroid(obj.Route(i,1)),obj.centroidXYZ(obj.Route(i,1),:));
                        geodesic_propagate(obj.Algorithm, source_points);
                        destination = geodesic_create_surface_point('face',obj.iCentroid(obj.Route(i,2)),obj.centroidXYZ(obj.Route(i,2),:));
                        obj.Path{i} = geodesic_trace_back(obj.Algorithm, destination);
                    end
            end
        end
        
        function getPathX(obj)
            % Get the co-ordinates of the computed geodesic path as they
            % cross the mesh edges
            delete(obj.PathHandles)
            x = [];
            y = [];
            z = [];
            for i = 1:size(obj.Route,1)
                [xTemp,yTemp,zTemp] = extract_coordinates_from_path(obj.Path{obj.Route(i,1)});
                x = [x; xTemp];
                y = [y; yTemp];
                z = [z; zTemp];
            end
            obj.PathX = [x y z];
        end
        
        function plotPath(obj)
            % Plot the geodesic path
            clear obj.PathHandles
            hold on
            obj.PathHandles = plot3( obj.PathX(:,1)*1.001 ...
                ,obj.PathX(:,2)*1.001 ...
                ,obj.PathX(:,3)*1.001 ...
                ,'k-s' ...
                , 'LineWidth', 2 ...
                , 'markersize', 10 ...
                );
        end
        
        function deleteGeodesic(obj)
            % Delete the geodesic mesh and algorithm from memory
            geodesic_delete()
        end
        
        function valueChangedFcn(obj, ~, ~)
            % If a DataTip location has changed, recompute the geodesic
            % path
            disp('in valueChangedFcn')
            obj.getDataTips();
            obj.labelDataTips();
            obj.getDataTipCoords();
            obj.initialiseGeodesicCalculation();
            obj.computePointOrder();
            if strcmpi(obj.Method, 'centroid')
                obj.convertDataTipCoordsToCentroids();
            end
            obj.computePath();
            obj.getPathX();
            obj.plotPath()
        end
        
        function getPathNodeDistances(obj)
            % Get closest vertex indices and distances to these vertices
            % for every location on the computed geodesic path
            [vertices, distances] = findclosestvertex(obj.Mesh, obj.PathX);
            obj.PathNodeDistances = [vertices distances];
        end
        
        function eCross = edgesIntersectingPath(obj)
            % Identify the mesh edges that the geodesic path crosses
            [~,~,edgeIndex,~] = triclosestedge(obj.PathX,obj.Mesh);
            e = edges(obj.Mesh);
            eCross = e(edgeIndex,:);
            eCross = unique(eCross,'rows');
            
            % Highlight the mesh edges crossed by the geodesic path
            hold on
            for i = 1:size(eCross,1)
                X = obj.Mesh.Points(eCross(i,:),:);
                plot3(X(:,1)*1.001,X(:,2)*1.001,X(:,3)*1.001,'r-','LineWidth',2);
            end
            set(obj.Surface, 'edgecolor', [.5 .5 .5])
        end
        
        function iFaces = identifyFacesIntersectingPath(obj)
           obj.eCross = obj.edgesIntersectingPath();
           iFaces = triedges2faces(obj.Mesh, obj.eCross);
           trNew = triangulation(obj.Mesh.ConnectivityList(iFaces,:), obj.Mesh.Points);
           trNew = repack(trNew);
           %trisurf(trNew);
        end
        
        function loadDataTipPoints(obj)
            load('datatiplocations.mat');
            for i = 1:size(dataTipX,1)
                datatip(obj.Surface, dataTipX(i,1), dataTipX(i,2), dataTipX(i,3));
            end
            obj.getDataTips();
        end
        
        function loadState(obj)
            obj.loadDataTipPoints
            load('holecutterstate.mat');
            obj.Method = Method;
            obj.XYZ = XYZ;
            obj.iV = iV;
            obj.centroidXYZ = centroidXYZ;
            obj.iCentroid =  iCentroid;
            obj.Route = Route;
            obj.Algorithm = Algorithm;
            obj.Path = Path;
            obj.PathX = PathX;
            obj.PathNodeDistances = PathNodeDistances;
            obj.plotPath();
        end
        
        function saveDataTipLocations(obj)
            for i = 1:numel(obj.DataTips)
                dataTipX(i,1) = obj.DataTips(i).X;
                dataTipX(i,2) = obj.DataTips(i).Y;
                dataTipX(i,3) = obj.DataTips(i).Z;
            end
            save('datatiplocations', 'dataTipX', '-v7.3');
        end
        
        function saveState(obj)
            obj.saveDataTipLocations();
            Method = obj.Method;
            XYZ = obj.XYZ;
            iV = obj.iV;
            centroidXYZ = obj.centroidXYZ;
            iCentroid = obj.iCentroid;
            Route = obj.Route;
            Algorithm = obj.Algorithm;
            Path = obj.Path;
            PathX = obj.PathX;
            PathNodeDistances = obj.PathNodeDistances;
            save('holecutterstate' ...
                , 'Method', 'XYZ', 'iV', 'centroidXYZ', 'iCentroid' ...
                , 'Route', 'Algorithm', 'Path', 'PathX', 'PathNodeDistances');
        end
        
        function erase(obj)
            % erase everything
            obj.getDataTips()
            try
                obj.deleteGeodesic();
            catch
                disp('HOLECUTTER/There is no geodesic library to unload');
            end
            delete(obj.DataTips)
            delete(obj.PathHandles)
            obj.DataTips = [];
            obj.PathHandles = [];
            obj.XYZ = [];
            obj.iV = [];
            obj.centroidXYZ = [];
            obj.iCentroid = [];
            obj.Route = [];
            obj.Algorithm = [];
            obj.Path = [];
            obj.PathX = [];
            obj.PathNodeDistances = [];
        end
        
        function computeCutOut(obj)
            obj.getDataTips();
            obj.getDataTipCoords();
            obj.initialiseGeodesicCalculation();
            obj.computePointOrder();
            obj.convertDataTipCoordsToCentroids();
            obj.Method = 'centroid';
            obj.computePath();
            obj.getPathX();
            obj.plotPath();
            obj.deleteGeodesic();
        end
        
        function newSurfaces = cutMesh(obj)
            % returns the newSurfaces as a cell array of triangulations.
            % The largest surface will always be newSurfaces{1}.

            iFaces = obj.identifyFacesIntersectingPath;
            trNew = triangulation(obj.Mesh.ConnectivityList(~iFaces,:), obj.Mesh.Points);
            
            % start with any edge, might as well be the first one, in eCross
            trTemp{1} = repack(getallattachedsurface(trNew, obj.eCross(1,1), 'vertex', 'verbose', false));
            trTemp{2} = repack(getallattachedsurface(trNew, obj.eCross(1,2), 'vertex', 'verbose', false));
            
            areas = [sum(triarea(trTemp{1})) sum(triarea(trTemp{2}))];
            [~,i] = max(areas);
            newSurfaces{1} = trTemp{i};
            newSurfaces{2} = trTemp{abs(i - 2) + 1};
            
            trisurf(newSurfaces{1}, 'facecolor', 'r');
        end
        
    end
end
