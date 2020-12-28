function areas = voltageHistogramAnalysis( userdata, varargin )
% VOLTAGEHISTOGRAMANALYSIS Performs voltage histogram analysis
%
% Usage:
%   areas = voltageHistogramAnalysis( userdata, varargin )
% Where:
%   userdata  - see importcarto_mem
%   areas     - he chamber areas within each of the voltage thresholds
%
% VOLTAGEHISTOGRAMANALYSIS accepts the following parameter-value pairs
%   'method'    {'map'}|'egm'
%   'type'      {'bip'}|'uni'
%   'threshold' n x 2 matrix of threshold values, default:
%             [ 0.01 0.11; 0.11 0.21; 0.21 0.30; 0.30 0.40; 0.40 0.50 ]
%   'plot'      {false} | true
%   'colors'      {colorBrewer colors r, y, g, b, p}|
%
% VOLTAGEHISTOGRAMANALYSIS displays a histogram of voltages coloured according
% to voltages, threshold. If 'method' is set to 'egm' then the bipolar 
% voltage is first interpolated from the bipolar electrogram data (stored 
% in userdata.electric). If 'type' is set to 'uni' then unipolar voltages 
% are used.
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

nStandardArgs = 1; % UPDATE VALUE
method = 'map';
type = 'bip';
threshold = [ 0.01 0.11; 0.11 0.21; 0.21 0.30; 0.30 0.40; 0.40 0.50 ];
colors = [colorBrewer('r'); colorBrewer('y'); colorBrewer('g'); colorBrewer('b'); colorBrewer('p')];
plot = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'method'
                method = varargin{i+1};
            case 'type'
                type = varargin{i+1};
            case 'threshold'
                threshold = varargin{i+1};
            case 'plot'
                plot = varargin{i+1};
            case 'colors'
                colors = varargin{i+1};
        end
    end
end
% TODO: check format of input values for each parameter are correct

% TODO: separate histogram analysis code so that it can be re-used for
% other purposes e.g. local activation time or conduction velocity
for i = 1:size(threshold, 1)
    [areas(i) voltages, ~, tr2{i}] =  getLowVoltageArea(userdata ...
        , 'threshold', threshold(i,:) ...
        , 'method', method ...
        , 'type', type ...
        );
end

if plot  
    % plot the histogram
    hold on;
    % plot bars for each of the thresholds specified
    for i = 1:size(threshold, 1)
        y1 = 0;
        y2 = areas(i);
        x1 = threshold(i,1);
        x2 = threshold(i,2);
        patch([x1 x1 x2 x2],[y1 y2 y2 y1], colors(i,:) ...
            , 'edgecolor', 'none' ...
            );
    end
    % calculate the number of other bars we need
    meanWidth = mean(range(threshold,2));
    remainingVoltages = voltages(voltages>threshold(end));
    nbins = floor(range(remainingVoltages) / meanWidth);
    for i = 1:nbins
        vstart = threshold(end)+(i-1)*meanWidth;
        vend = threshold(end)+i*meanWidth;
        area =  getLowVoltageArea(userdata ...
        , 'threshold', [vstart vend] ...
        , 'method', method ...
        , 'type', type ...
        );
        y1 = 0;
        y2 = area;
        x1 = vstart;
        x2 = vend;
        patch([x1 x1 x2 x2],[y1 y2 y2 y1], [.5 .5 .5] ...
            , 'edgecolor', 'none' ...
            );
    end 
    xlabel('Voltage (mV)');
    ylabel('Area (cm^2)');
    set(gcf, 'color', 'w')
    set(gca, 'fontsize', 14); 
    
        % plot the surface
    figure
    drawMap(userdata, 'type', 'none', 'orientation', 'pa');
    hold on
    for i = 1:size(threshold, 1)
        if ~isempty(tr2{i})
            trisurf(tr2{i} ...
                , 'facecolor', colors(i,:) ...
                , 'edgecolor', 'none' ...
                , 'SpecularStrength', 0 ...
                );
        end
    end
end
