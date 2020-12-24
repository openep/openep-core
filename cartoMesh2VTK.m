function tr = cartoMesh2VTK()
% CARTOMESH2VTK Converts a Carto mesh file to VTK file
%
% Usage:
%   tr = cartoMesh2VTK('openfile')
% Where:
%   tr,  - a TriRep object
%
% CARTOMESH2VTK Converts a Carto3 mesh to a VTK file and returns a TriRep
% object
%
% Author: Steven Williams (2015) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
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

% output the mesh to a temporary VTK
[fileName, pathName] = uigetfile('*.mesh', 'Select a Carto3 .mesh file');
tr = read_meshfile([pathName filesep fileName]);
fileName2 = [fileName(1:end-4) 'vtk'];
writeTriRep2VTK(repack(tr), [], 'outputfile', [pathName filesep fileName2]);

end
