function plotAblationArea(userdata)
% PLOTABLATIONAREA Adds the ablation area to the current figure
%
% Usage:
%   plotAblationArea(userdata)
% Where:
%   userdata - see importcarto_mem.m
%
% PLOTABLATIONAREA Requires a userdata structure which contains `.rfindex` as
% its input, which can be created using `importvisitag.m`
%
% Author: Steven Williams (2020) (Copyright)
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
%
% See also, importvisitag.m, getAblationArea.m, plotVisitags.m

[~, ~, trAbl] = getAblationArea(userdata);
trisurf(trAbl, 'facecolor', 'y', 'edgecolor', 'none');

end