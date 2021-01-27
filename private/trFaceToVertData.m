function [ vertexData ] = trFaceToVertData( tr, faceData )
% TRFACETOVERTDATA converts face to vertex data for a TriRep object.
% Usage:
%   vertexData = trFaceToVertData(tr, faceData)
% Where:
%   tr - is a TriRep object with the following dimensions:
%           tr.Triangulation    q * 3
%           tr.X                p * 1
%   faceData - is an p*1 array of scalar data (by face)
%   vertexData - is a q*1 array of scalar data (by vertex)
%
% TRFACETOVERTDATA detailed description.
%
% Author: Steven Williams (2012)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% Use warnings with the following format
%    disp('TRFACETOVERTDATA: warning.')
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

[~, locs] = tricentroid(tr);

F = scatteredInterpolant(locs(:,1), locs(:,2), locs(:,3) ...
    , faceData ...
    , 'linear' ... % interpolation {nearest|linear|natural}
    , 'linear' ... % extrapolation {nearest|linear|none}
    );

vertexData = F(tr.X);