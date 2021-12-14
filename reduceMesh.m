function newUserdata = reduceMesh( userdata, R )
% REDUCEMESH reduces the number of triangles in the OpenEP dataset surface
%
% Usage:
%   newUserdata = reduceMesh( userdata )
% Where:
%   userdata - the input OpenEP dataset, see https://openep.io/data/
%   R - see, reducepatch.m
%
% REDUCEMESH reduces the number of faces.
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

F = getFaces(userdata);
V = getVertices(userdata);

[nf, nv] = reducepatch(F, V, R);

newUserdata.surface.triRep = TriRep(nf, nv);

end