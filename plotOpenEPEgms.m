function hFig = plotOpenEPEgms( userdata, varargin )
% PLOTOPENEPEGMS Plot eletrograms from OpenEP data
%
% Usage:
%   [ hFig ] = plotOpenEPEgms( userdata, varargin )
%
% Where:
%   userdata  - see importcarto_mem
%   hFig - a handle to the plotted figure
%
% PLOTOPENEPEGMS accepts the following parameter-value pairs
%   'iEgm'     {:}|[a:b]
%           an array indexing into userdata.electric.egm such that
%          %   'iEgm'     {:}|[a:b]
%           an array indexing into userdata.electric.egm such that
%           userdata.electric.egm(iEgmArray,:) are selected for plotting
%           To convert from Carto point numbers to iEgmArray use
%           getIndexFromCartoPointNumber.
%   'range'         {'window'}|'all'
%           By default ('window') only the electrogram within the window of
%           interest is drawn (±buffer). By specifying 'all' the entire point
%           electrogram is drawn.
%   'buffer'        {50}|double
%           The time before and after the window of interest to draw. By
%           default, 20ms, but can be changed by setting 'buffer' to an
%           alternative value.
%   'egmtype'   'bip'|'uni'|{'bip-uni'}
%           Whether to plot only the bipolar electrograms, only the
%           unipolar electrograms or both
%   'reference' 'off'|{'on'}
%           Whether to plot the reference channel, off by default userdata.electric.egm(iEgmArray,:) are selected for plotting
%           To convert from Carto point numbers to iEgmArray use
%           getIndexFromCartoPointNumber.
%
% PLOTOPENEPEGMS is a wrapper function for plotElectrograms.
%
% Author: Steven Williams (2017) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%   SW 2020: Modified for incorporation into OpenEP
%
% See also PLOTELECTROGRAMS, GETINDEXFROMCARTOPOINTNUMBER,
% GETELECTROGRAMSATPOINTS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% plotOpenEPEgms(userdata, 'iegm', getIndexFromCartoPointNumber(userdata,1))
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 2; % UPDATE VALUE
iEgmArray = ':';
range = 'window';
egmtype = 'bip-uni';
buffer = 50;
reference = 'on';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'iegm'
                iEgmArray = varargin{i+1};
            case 'axis'
                hAx = varargin{i+1};
            case 'window'
                range = varargin{i+1};
            case 'buffer'
                buffer = varargin{i+1};
            case 'egmtype'
                egmtype = varargin{i+1};
            case 'reference'
                reference = varargin{i+1};
        end
    end
end
%TODO: Check validity of input measurements

% Create the cell array of electrograms
if strcmpi(iEgmArray, ':')
    iEgmArray = 1:numel(userdata.electric.tags);
end
nEgms = numel(iEgmArray);
for iEgm = 1:nEgms
    switch egmtype
        case 'bip'
            egmTraces{iEgm} = userdata.electric.egm(iEgmArray(iEgm),:);
            egmNames{iEgm} = char(userdata.electric.names(iEgmArray(iEgm),:));
            acttime(iEgm) = userdata.electric.annotations.mapAnnot(iEgmArray(iEgm));
            egmColors{iEgm} = 'b';
        case 'uni'
            for i = 1:2
                egmTraces{i,iEgm} = userdata.electric.egmUni(iEgmArray(iEgm),:,i);  
                acttime(i,iEgm) = userdata.electric.annotations.mapAnnot(iEgmArray(iEgm));
                egmNames{i,iEgm} = char(userdata.electric.names(iEgmArray(iEgm),:));
                egmColors{i,iEgm} = 'g';
            end
        case 'bip-uni'
            i = 1;
            egmTraces{i,iEgm} = userdata.electric.egm(iEgmArray(iEgm),:);
            egmNames{i,iEgm} = char(userdata.electric.names(iEgmArray(iEgm),:));
            acttime(i,iEgm) = userdata.electric.annotations.mapAnnot(iEgmArray(iEgm));
            egmColors{i,iEgm} = 'b';
            for i = 1:2
                egmTraces{i+1,iEgm} = userdata.electric.egmUni(iEgmArray(iEgm),:,i);  
                acttime(i+1,iEgm) = userdata.electric.annotations.mapAnnot(iEgmArray(iEgm));
                egmNames{i+1,iEgm} = char(userdata.electric.names(iEgmArray(iEgm),:));
                egmColors{i+1,iEgm} = 'g';
            end
    end
    if strcmpi(reference, 'on')
        egmTraces{end+1,iEgm} = userdata.electric.egmRef(iEgmArray(iEgm),:) / 20;
        acttime(end+1,iEgm) = userdata.electric.annotations.referenceAnnot(iEgm);
        egmNames{end+1,iEgm} = 'Ref';
        egmColors{end+1,iEgm} = 'r';
    end
end
egmTraces = flipud(egmTraces(:));
acttime = flipud(acttime(:));
egmNames = flipud(egmNames(:));

egmColors = flipud(egmColors(:));

% Work out the range
switch range
    case 'all'
        range = [NaN NaN];
    case 'window'
        range = userdata.electric.annotations.woi(iEgmArray(1),:) ...
            + userdata.electric.annotations.referenceAnnot(iEgmArray(1));
        buffer = [-buffer buffer];
        range = range + buffer;
end

% Finally, use plotElectrograms to actually draw the electrograms
hFig = plotElectrograms(egmTraces ...
    , 'egmNames', egmNames ...
    , 'range', range ...
    , 'acttime', acttime ...
    , 'egmColors', egmColors ...
    , 'separation', 7.5 ...
    , 'autogain', true ...
    , 'title', 'off' ...
    );


end