function [ hFig ] = plotprecision_modelgroups( tr, dxgeo )
% PLOTPRECISION_MODELGROUPS draws anatomies
%
% Usage:
%   [ hFig ] = plotprecision_modelgroups( tr )
% Where:
%   tr  - a structure containing triangulation objects, see processprecision_modelgroups.m
%   dxgeo - see loadprecision_modelgroups.m
%   hFig - a structure of figure handles
%
% PLOTPRECISION_MODELGROUPS Detailed description goes here
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

for i = 1:numel(tr)
    hFig{i} = figure;
    trisurf(tr{i}, 'facecolor', 'g', 'edgecolor', 'none');
    title(dxgeo(i).name);
    axis equal vis3d
    cameraLight;
end

end