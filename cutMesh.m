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
hc = HoleCutter(getMesh(userdata, 'type', 'triangulation', 'repack', true), hT);

disp('hc.computeCutOut')
hc.computeCutOut

disp('hc.cutMesh')
newSurfaces = hc.cutMesh()

% now divide the userdata into two different sections
tf = local_splitSurfaceData(getVertices(userdata), newSurfaces{1}.Points);

userdata1 = userdata;
userdata1.surface.triRep = newSurfaces{1};
userdata1.surface.isVertexAtRim = local_isVertexAtRim(newSurfaces{1});
userdata1.surface.act_bip = userdata.surface.act_bip(tf,:);
userdata1.surface.uni_imp_frc = userdata.surface.uni_imp_frc(tf,:);

tf = local_splitSurfaceData(getVertices(userdata), newSurfaces{2}.Points);

userdata2 = userdata;
userdata2.surface.triRep = newSurfaces{2};
userdata2.surface.isVertexAtRim = local_isVertexAtRim(newSurfaces{2});
userdata2.surface.act_bip = userdata.surface.act_bip(tf,:);
userdata2.surface.uni_imp_frc = userdata.surface.uni_imp_frc(tf,:);

newUserdata{1} = userdata1;
newUserdata{2} = userdata2;

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

    function  tf = local_splitSurfaceData(pointSetOriginal, pointSetNew)
        % find the set of indices into pointSetOriginal for the vertices
        % which are closest to (in this case the same as) pointSetNew.

        tf = ismember(pointSetOriginal, pointSetNew, 'rows');
    end


end