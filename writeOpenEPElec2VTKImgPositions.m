function writeOpenEPElec2VTKImgPositions(Y, colourIndex, colourMap)
% WRITEOPENEPELEC2VTKIMGPOSITIONS Writes a VTK file containing eletrode
% positions mapped to imaging data
%
% Usage:
%   writeOpenEPElec2VTKImgPositions(Y, colourIndex, colourMap, filename)
% Where:
%   Y             - 3D Cartesian co-ordinates of the electrode positions in imaging space
%
% WRITEOPENEPELEC2VTKIMGPOSITIONS Can be called to write a VTK file after
% Y, colourIndex and colourMap have been created by TRANSOPENEP2IMGPOS. The
% colour values, which are based on a 2D colour scale with values selected
% by the UACs, are written directly into the VTK file. To visualise these
% points in Paraview, select "PointGuassian" and untick "Map scalars" under
% the "Scalar Coloring" section. (Note, may need to turn on advance options
% by clicking the gear icon in the Paraview property editor).
%
% Author: Steven Williams (2021)
% Modifications -
%
% See also: TRANSOPENEP2IMGPOS
%
% Info on Code Testing:
% ---------------------------------------------------------------
%  
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

hVtk = VTKWriter('openfile');
hVtk.Points = Y'/1000;

C = colourMap(colourIndex,:);
C = C/256;
C(:,4) = 1;
hVtk.ColorScalars = C';

x = Y(:,1);
y = Y(:,2);
z = Y(:,3);
DT = delaunayTriangulation(x,y,z);
[K,~] = convexHull(DT);

hVtk.Cells = K';

hVtk.writeVTK;