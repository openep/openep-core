function cellLabels = labelOpenEPMapFromImgMap( userdata, mOpenEP3D )
% LABELOPENEPMAPFROMIMGMAP Label OpenEP data from an imaging mesh with
% different resolution
%
% Usage:
%   cellLabels = labelOpenEPMapFromImgMap( userdata, mOpenEP3D )
% Where:
%   userdata  - see importcarto_mem
%   mOpenEP3D - a mesh structure with .Pts and .Tri fields, where the
%               fourth column of .Tri is the labels
%   cellLabels - the labels, referencing into userdata.surface
%
% LABELOPENEPMAPFROMIMGMAP works by finding the closest cell in mOpenEP3D
% to each cell in userdata, and assigning the label from mOpenEP3D to the
% cell in userdata. 
%
% TODO: Deal with vertex labelling
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% 
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

trOpenEP = userdata.surface.triRep;
[~, trCentroid_OpenEP] = tricentroid(trOpenEP); 

trImg = getTriangulationFromMeshStruct( mOpenEP3D, 'region', ':', 'scale', 'um', 'repack', 'false');
[~, trCentroid_Img] = tricentroid(trImg);

iCentroid = findclosestvertex(trCentroid_Img, trCentroid_OpenEP);

cellLabels = mOpenEP3D.Tri(iCentroid,4);

end