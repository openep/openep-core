function [ userdata1, userdata2 ] = cutMesh( userdata )
% CUTMESH Cut the mesh into two parts around any given continuous path
%
% Usage:
%   [ userdata1, userdata2 ] = cutMesh( userdata )
% Where:
%   userdata - the input OpenEP dataset, see https://openep.io/data/
%   userdata1 - the first output
%   userdata2 - the second output
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
newSurface = hc.cutMesh()

userdata1 = userdata;
userdata2 = userdata;

end