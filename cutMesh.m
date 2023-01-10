function [ newUserdata, newSurfaces ] = cutMesh( userdata )
% CUTMESH Cut the mesh into two parts around any given continuous path
%
% Usage:
%   [ userdata1, userdata2 ] = cutMesh( userdata )
% Where:
%   userdata - the input OpenEP dataset, see https://openep.io/data/
%   newUserdata
%
% CUTMESH can, for example, be used to create a valve cut-out, crop the
% aortic geometry, remove pulmonary veins. CUTMESH uses a HoleCutter
% object to perform mesh manipulation. CUTMESH also depends on geodesic
% path calculations.
%
% Author: Steven Williams (2021)
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

hT = drawMap(userdata, 'type', 'none');
set(hT, 'facealpha', 0.8);

% set the tool tip format
hDt = datatip(hT, 'dataindex', 1);
hT.DataTipTemplate.DataTipRows(2) = [];
hT.DataTipTemplate.DataTipRows(2) = [];
hT.DataTipTemplate.DataTipRows(1).Label = 'x';
hT.DataTipTemplate.DataTipRows(1).Value = [];
delete(hDt);

% prompt user to select tool tips
disp('Select points on the mesh and press any key when done')
pause;

% perform mesh cutting
hc = HoleCutter(getMesh(userdata, 'type', 'triangulation', 'limitToTriangulation', true), hT);

disp('hc.computeCutOut')
hc.computeCutOut

disp('hc.cutMesh')
newSurfaces = hc.cutMesh()

% divide the userdata into two different sections, dealing with surface data
tf = local_labelSurfaceData(getVertices(userdata), newSurfaces{1}.Points);
userdata1 = local_splitSurfaceData(userdata, newSurfaces{1}, tf);

tf = local_labelSurfaceData(getVertices(userdata), newSurfaces{2}.Points);
userdata2 = local_splitSurfaceData(userdata, newSurfaces{2}, tf);

% repack the surface data
newUserdata{1} = repackUserdata(userdata1);
newUserdata{2} = repackUserdata(userdata2);

% deal with electric data
electricIds = local_labelElectricData(userdata, newUserdata);
newUserdata{1} = local_splitElectricData(newUserdata{1}, electricIds, 2);
newUserdata{2} = local_splitElectricData(newUserdata{2}, electricIds, 1);

    function isVertexAtEdge = local_isVertexAtRim(tr)
        if isa(tr, 'TriRep')
            xyz = tr.X;
        elseif isa(tr, 'triangulation')
            xyz = tr.Points;
        else
            error('OPENEP/CUTMESH: Wrong data format - must be TriRep or triangulation');
        end
        isVertexAtEdge = false(size(xyz,1),1);
        vFree = freeBoundary(tr);
        vFree = vFree(:);
        [~, iFB, ~] = unique(vFree);
        vFree(iFB) = [];    % we have wiped out all vertices which are only listed once
        isVertexAtEdge(vFree) = true;
    end

    function  tf = local_labelSurfaceData(pointSetOriginal, pointSetNew)
        % find the set of indices into pointSetOriginal for the vertices
        % which are closest to (in this case the same as) pointSetNew.
        tf = ismember(pointSetOriginal, pointSetNew, 'rows');
    end

    function electricIds = local_labelElectricData(userdata, newUserdata)
        % index the electric data by finding the surface closest to
        % each electric point.
        [~, distances1] = findclosestvertex(getVertices(newUserdata{1}), getElectrogramX(userdata));
        [~, distances2] = findclosestvertex(getVertices(newUserdata{2}), getElectrogramX(userdata));
        allDistances = [distances1 distances2];
        [~,electricIds] = min(allDistances, [], 2);
    end

    function userdata1 = local_splitSurfaceData(userdata, newSurface, tf)
        userdata1 = userdata;
        userdata1.surface.triRep = newSurface;
        userdata1.surface.isVertexAtRim = local_isVertexAtRim(newSurface);
        if ~isempty(userdata.surface.act_bip)
            userdata1.surface.act_bip = userdata.surface.act_bip(tf,:);
        end
        if ~isempty(userdata.surface.uni_imp_frc)
            userdata1.surface.uni_imp_frc = userdata.surface.uni_imp_frc(tf,:);
        end
    end

    function userdata = local_splitElectricData(userdata,electricIds,id)
        % subdivide the electric data by finding the surface closest to
        % each electric point.
        userdata.electric.tags(electricIds==id,:) = [];
        userdata.electric.names(electricIds==id,:) = [];
        userdata.electric.electrodeNames_bip(electricIds==id,:) = [];
        userdata.electric.egmX(electricIds==id,:) = [];
        userdata.electric.electrodeNames_uni(electricIds==id,:) = [];
        userdata.electric.egmUniX(electricIds==id,:,:) = [];
        userdata.electric.egmUni(electricIds==id,:,:) = [];
        userdata.electric.egmRef(electricIds==id,:) = [];
        userdata.electric.ecg(electricIds==id,:) = [];
        userdata.electric.annotations.woi(electricIds==id,:) = [];
        userdata.electric.annotations.referenceAnnot(electricIds==id,:) = [];
        userdata.electric.annotations.mapAnnot(electricIds==id,:) = [];
        userdata.electric.voltages.bipolar(electricIds==id,:) = [];
        userdata.electric.voltages.unipolar(electricIds==id,:) = [];
        userdata.electric.impedances.time(:,electricIds==id) = [];
        userdata.electric.impedances.value(:,electricIds==id) = [];
        userdata.electric.egmSurfX(electricIds==id,:) = [];
        userdata.electric.barDirection(electricIds==id,:) = [];
    end
end