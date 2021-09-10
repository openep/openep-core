function [normals, userdata] = getNormals( userdata )
% GETNORMALS computes or returns surface normals on OpenEP anatomy
%
% Usage:
%   normals = computeNormals( userdata )
%   [normals, userdata] = computeNormals( userdata )
% Where:
%   userdata  - see importcarto_mem
%   normals - the surface normals at each vertex
%
% COMPUTENORMALS Returns the normals but also returns a new userdata
% structure with normals saved
%
% Author: Steven Williams (2021) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% faces = getFaces( userdata )
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

if isfield(userdata.surface, 'normals')
    if isempty(userdata.surface.normals)
        V = getVertices(userdata);
        T = getFaces(userdata);
        normals = local_computeNormals(V,T);
    else
        normals = userdata.surface.normals;
    end
else
    V = getVertices(userdata);
    T = getFaces(userdata);
    normals = local_computeNormals(V,T);
end

userdata.surface.normals = normals;

    function N = local_computeNormals(V, T)
        N = compute_vertex_normals(V, T, 6, 'norm');
    end

end
