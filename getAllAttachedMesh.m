function newUserdata = getAllAttachedMesh( userdata )
% GETALLATTACHEDMESH Ensures the mesh in userdata is continuous
%
% Usage:
%   newUserdata = getAllAttachedMesh( userdata )
% Where:
%   userdata - the input OpenEP dataset, see https://openep.io/data/
%   newUserdata - the output OpenEP dataset
%
% GETALLATTACHEDMESH ensures that the surface in userdata is continuous and
% returns a new userdata with the surface modified.
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
disp('Select a single point on the mesh and press any key when done')
pause;

% Get the co-ordinate of the point and the vertex index
dataTip = get(hT, 'children');
XYZ(1,1:3) = [dataTip(1).X dataTip(1).Y dataTip(1).Z];
trSurface = getMesh(userdata, 'type', 'triangulation', 'repack', true)
iV(1) = findclosestvertex(trSurface, XYZ(1,:));

% Get all the attached surface
newTriRep = getallattachedsurface(trSurface, iV, 'vertex');

% Store the new surface
newUserdata = userdata;
newUserdata.surface.triRep = newTriRep;

end