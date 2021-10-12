function [ tr ] = processprecision_modelgroups( dxgeo )
% PROCESSPRECISION_MODELGROUPS Creates triangulations from dxgeo
%
% Usage:
%   [ tr ] = processprecision_modelgroups( dxgeo )
% Where:
%   dxgeo  - see loadprecision_modelgroups.m
%   tr  - a structure containing triangulation objects
%
% PROCESSPRECISION_MODELGROUPS Detailed description goes here
%
% Author: Steven Williams (2017)
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

for i = 1:numel(dxgeo)
    tr{i} = triangulation(dxgeo(i).triangles, dxgeo(i).vertices);
end

end