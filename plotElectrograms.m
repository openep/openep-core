function hFig = plotElectrograms(egmTraces, varargin)
% PLOTELECTROGRAMS Draws any number of electrograms in line.
%
% Usage:
%   hFig = PLOTELECTROGRAMS(egmTraces)
%   hFig = PLOTELECTROGRAMS(egmTraces, varargin)
%
% PLOTELECTROGRAMS parameter-value pairs can be passed in as follows:
% 'egmNames', cell array of electrogram names
% 'range', two value vector [xmin xmax]
% 'sampleRate', the sample rate in Hz
% 'paperSpeed', the paper speed for plotting
% 'separation', the separation between electrograms, default is 5mV
% 'clipping', two value vector [clipmin clipmax] in mV
% 'gain', the gain, default is 2
% 'autogain', if true the electrograms are automatically scaled, default false
% 'acttime', a vector of activation times to be plotted as red crosses
% 'egmColors', cell array of colors. Same size as egmNames
%
% INSTRUCTIONS FOR A4 EGM PAGE
% 1. Save figure as an .EPS file.
% 2. Open in Corel draw and export as a .PDF file.
%
% INSTRUCTIONS FOR TIMING SINGLE BEAT TRACES
% The default position gives a width of 300pixels, which corresponds to an
% x axis width of 68.44mm. The number of samples to plot in this
% width (ie width of 'range' parameter) is given by:
%   n = (sampleRate / paperSpeed) * 68.44
% For example, at 1000Hz range should be specified as having width:
%   25mm/s  - 2738  [-1119 1619]
%   100mm/s - 684
%   200mm/s - 342
%
% Author: Steven Williams (2012) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%   SW 2020: Modified function to incorporate into OpenEP.
%
% See also PLOTOPENEPEGMS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% For an A4 page:
%     pageSize = [5 0 29.7 21.0]; %A4 landscape
%     plotElectrogram(egm ...
%         , 'range', axisRange ...
%         , 'paperSpeed', hBg.paperSpeed ...
%         , 'position', pageSize ...
%         , 'units', 'centimeters' ...
%         , 'name', 'Exported from BardGUI' ...
%         , 'sampleRate', hBg.hBard.SampleRate ...
%         , 'autogain', true ...
%         );
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% Validate input
p = inputParser;
p.addRequired('egmTraces', @iscell);
p.addParamValue('egmNames', {}, @iscellstr);
p.addParamValue('range', [NaN NaN], @(x)isnumeric(x) && numel(x)==2);
p.addParamValue('sampleRate', NaN, @(x)isnumeric(x) && numel(x)==1);
p.addParamValue('paperSpeed', NaN, @isnumeric);
p.addParamValue('units', 'pixels', @(x)ischar(x) && strcmpi(x, 'pixels') || strcmpi(x, 'centimeters'));
p.addParamValue('position', [200 200 300 600], @(x)isnumeric(x) && numel(x)==4);
p.addParamValue('name', 'name', @isstr);
p.addParamValue('separation', 5, @(x)isnumeric(x) && numel(x)==1);
p.addParamValue('clipping', [-4 4], @(x)isnumeric(x) && numel(x)==2);
p.addParamValue('gain', 2, @(x)isnumeric(x) && numel(x)==1);
p.addParamValue('autogain', false, @islogical);
p.addParamValue('center', false, @islogical);
p.addParamValue('acttime', 0, @isnumeric);
p.addParamValue('egmColors', {}, @iscellstr);
p.addParamValue('axis', []);
p.addParamValue('title', 'on', @isstr);
p.parse(egmTraces, varargin{:});
inputs = p.Results;

egmTraces = inputs.egmTraces;

numEgm = numel(egmTraces);
len = NaN(numEgm,1);
for i = 1:numEgm
    len(i) = length(egmTraces{i});
end
lengthEgm = nanmax(len);

paperSpeed = inputs.paperSpeed;
units = inputs.units;
position = inputs.position;
name = inputs.name;
separation = inputs.separation; %mV
clipping = inputs.clipping; %mV
gain = inputs.gain;
autogain = inputs.autogain;
center = inputs.center;
acttime = inputs.acttime;
ax = inputs.axis;
titleText = inputs.title;

if isempty(inputs.egmNames)
    egmNames = num2cell(1:numEgm);
else
    egmNames = inputs.egmNames;
end
if isnan(inputs.range)
    range = [1 lengthEgm];
else
    range = inputs.range;
end
if isnan(inputs.sampleRate)
    sampleRate = 1000;
    sampleRateTitle = 'NaN';
    xAxisLabel = 'samples';
else
    sampleRate = inputs.sampleRate;
    sampleRateTitle = num2str(sampleRate);
    xAxisLabel = 'time (s)';
end
if ~isempty(inputs.egmColors)
    egmColors = fliplr(inputs.egmColors);
else
    egmColors = cell(1,numel(egmNames));
    egmColors(:) = {'k'};
end

% Set up the figure
if isempty(ax)
    hFig = figure;
    set(hFig, 'color', 'w' ...
        , 'units', units ...
        , 'position', position ...
        , 'name', name ...
        );
end

% If necessary add zeros to the start or end of the trace
deltaTime = range(2) - range(1) + 1;
egm = NaN(floor(deltaTime), numEgm);
for i = 1:numEgm
    try
        egm(:,i) = egmTraces{i};
    catch MeX
        if center
            lengthDif = length(egm(:,i)) - length(egmTraces{i});
            egm(floor(lengthDif/2+1):floor(lengthDif/2+length(egmTraces{i})),i) = egmTraces{i};
        else
            if range(1) < 0 % fill the left of the trace with zeros
                egm(end-length(egmTraces{i})+1:end,i) = egmTraces{i};
                disp(['PLOTELECTROGRAM: ' MeX.message ' Adding zeros to electrogram: ' egmNames{i}]);
            else % fill the right of the trace with zeros
                egm(1:length(egmTraces{i}),i) = egmTraces{i};
                disp(['PLOTELECTROGRAM: ' MeX.message ' Adding zeros to electrogram: ' egmNames{i}]);
            end
        end
    end
end

% Apply gain, clipping and spacing
if (autogain)
    for i = 1:numEgm
       egmDelta = nanmax(egm(:,i)) - nanmin(egm(:,i));
       multiplier = 0.9*separation/egmDelta;
       egm(:,i) = egm(:,i) * multiplier;
    end
else
    egm = egm .* gain;
end
egm(egm>clipping(2)) = clipping(2);
egm(egm<clipping(1)) = clipping(1);
for i = 1:numEgm
    egm(:,i) = egm(:,i) + separation * i;
end

% Draw the electrograms
time = range(1):range(2);
time = time';
xValues = repmat(time, [1, numEgm]);
for i = 1:size(xValues,2)
    line(xValues(:,i), egm(time,i) ...
        , 'color', colorBrewer(egmColors{i}) ...
        , 'linewidth', 1 ...
        );
    
    yTickPos = separation:separation:separation*numEgm;
    yTickLabel = egmNames;
    
    set(gca, 'YTick', yTickPos ...
        , 'YTickLabel', yTickLabel ...
        , 'YLim', [0 separation*(numEgm+1)] ...
        , 'XLim', [range(1) range(2)] ...
        , 'Color', 'w' ...
        );
    %    , 'Position', [0.09 0.11 0.861 0.815] ... % Fine adjustment for speed
    
    
    % Draw the acttimes
    if ~acttime==0
        hold on
        yVal = separation:separation:numel(acttime)*separation;
        plot(acttime, yVal, 'r.', 'markersize', 20);
    end
    
    % Sort out the x axis labels
    xTick = get(gca, 'XTick');
    xTickLabel = xTick / (sampleRate/1000);
    set(gca, 'XTickLabel', round(xTickLabel));
    
    % Set the title and the x-axis label
    if strcmpi(titleText, 'off')
       % do nothing
    else
        title(['Paper Speed: ' num2str(paperSpeed) ' mm/s | Sample Rate: ' sampleRateTitle ' Hz']);
    end
    if isnan(sampleRate)
        xlabel('samples');
    else
        xlabel(xAxisLabel);
    end
    
    % Turn on the grid lines
    %set(gca, 'XGrid', 'on');
    
end