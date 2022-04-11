function plotVisitags(userdata, varargin)
% PLOTVISITAGS Displays ablation data for a case
%
% Usage:
%   plotVisitags(userdata)
% Where:
%   userdata - see importcarto_mem.m
%
% PLOTVISITAG accepts the following parameter-value pairs
%   'plot'     {'tags'} | 'grid' | 'both'
%       - specifies whether to show the tags, the grid, or both
%   'shell'    {'on'} | 'off'
%       - specifies whether to show the chamber shell
%   'colour'    {'r'} | colorspec | array     
%       - can be a string or colorspec specifying the color of all the
%       spheres
%       - can be an array of double values which is rendered as a
%       colorscale
%   'orientation'   
%       - see `drawMap.m`
%   'tagindices'    {':'} | array
%       - the indices of the visitags to display
%
% PLOTVISITAG Requires a userdata structure which contains .rfindex as
% its input, which can be created using importvisitag.m
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% plotVisitags(userdata)
% plotVisitags(userdata, 'plot', 'both', 'shell', 'off', 'orientation', 'ap')
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------
%
% See also, importvisitag.m, getAblationArea.m, plotAblationArea.m

nStandardArgs = 1; % UPDATE VALUE
toplot = 'tags';
shell = 'on';
colour = 'r';
orientation = 'pa';
indices = ':';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'plot'
                toplot = varargin{i+1};
            case 'shell'
                shell = varargin{i+1};
            case 'colour'
                colour = varargin{i+1};
            case 'orientation'
                orientation = varargin{i+1};
            case 'tagindices'
                indices = varargin{i+1};
        end
    end
end
if strcmpi(indices, ':')
    indices = 1:size(userdata.rfindex.tag.X,1);
end
%TODO add error checking for the input param/value pairs

% Create a figure
figure
hold on

% Plot the shell, if needed or turn the light on
if strcmpi(shell, 'on')
    drawMap(userdata, 'type', 'none', 'orientation', orientation);
else
    cameraLight();
end

% Plot the tags
switch toplot
    case 'tags'
        facealpha = 1;
    case 'grid'
        facealpha = 1;
    case 'both'
        facealpha = 0.5;
end

% convert colour values into colormap indices
cMap = color_shades({'white' 'mistyrose' 'red'});
if isnumeric(colour)
    cInd = round((colour-min(colour)+1)/range(colour) * (size(cMap,1)-1));
end
if strcmpi(toplot, 'tags') || strcmpi(toplot, 'both')

    for i = 1:numel(indices) %size(userdata.rfindex.tag.X,1)
        % get the colour
        if ischar(colour)
            C = colour;
        else
            C = cMap(cInd(i),:);
        end
        %X = userdata.rfindex.tag.X(i,:);
        X = userdata.rfindex.tag.X(indices(i),:);
        h = plotsphere(X(1), X(2), X(3), C, 3, 16);
        set(h ...
            , 'facealpha', facealpha ...
            , 'edgecolor', 'none' ...
            );
    end
end

% Plot the grid
if strcmpi(toplot, 'grid') || strcmpi(toplot, 'both')
    for i = 1:numel(userdata.rfindex.grid)
        for j = 1:numel(userdata.rfindex.grid{i})
            X = userdata.rfindex.grid{i}(j).X;
            plot3(X(1), X(2), X(3), '.' ...
                , 'color', [0 0 0] ...
                , 'markersize', 8 ...
                );
        end
    end
end

end
