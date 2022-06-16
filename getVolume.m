function [ volume ] = getVolume( userdata )
% GETVOLUME Calculates the volume of the chamber described in userdata
%
% Usage:
%   volume = getVolume(userdata)
% Where:
%   userdata  - see importcarto_mem
%   volume  - the volume, in cm^3
%
% GETVOLUME For details of the calculation see:
%   https://stackoverflow.com/questions/1406029/how-to-calculate-the-volume-of-a-3d-mesh-object-the-surface-of-which-is-made-up
%   http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
%   EFFICIENT FEATURE EXTRACTION FOR 2D/3D OBJECTS IN MESH REPRESENTATION
%
% Author: Steven Williams (2017) (Copyright)
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

Tri = getMesh(userdata).Triangulation;
X = getMesh(userdata).X;

volume = 0;

for i = 1:length(Tri)
    v1 = Tri(i,1);
    v2 = Tri(i,2);
    v3 = Tri(i,3);
    
    p1 = X(v1,:);
    p2 = X(v2,:);
    p3 = X(v3,:);
    
    volume = volume + local_signedVolumeOfTriangle(p1, p2, p3);
end

volume = volume / 1000;

    function vol = local_signedVolumeOfTriangle(p1, p2, p3)
        % p1, p2, p3 are the three vertices of a triangle
        v321 = p3(1) * p2(2) * p1(3);
        v231 = p2(1) * p3(2) * p1(3);
        v312 = p3(1) * p1(2) * p2(3);
        v132 = p1(1) * p3(2) * p2(3);
        v213 = p2(1) * p1(2) * p3(3);
        v123 = p1(1) * p2(2) * p3(3);
        
        vol = (1/6) * (-v321 + v231 + v312 - v132 - v213 + v123);
    end

end