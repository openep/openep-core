function [egmTraces, acttime, egmNames] = getEgmsAtPoints( userdata, varargin )
% GETEGMSATPOINTS Access eletrograms from stored in the OpenEP data format
%
% Usage:
%   [egmTraces, acttime, egmNames] = getEgmsAtPoints( userdata, varargin )
%
% Where:
%   userdata  - see importcarto_mem
%   egmTraces - cell array of the requested electrograms
%   acttime   - cell array of activation times 
%   egmNames  - names of the electrograms
%
% GETEGMSATPOINTS accepts the following parameter-value pairs
%   'iEgm'     {:}|[a:b]
%           an array indexing into userdata.electric.egm such that
%           userdata.electric.egm(iEgmArray,:) are selected for plotting
%           To convert from Carto point numbers to iEgmArray use
%           getIndexFromCartoPointNumber.
%   'egmtype'   'bip'|'uni'|{'bip-uni'}
%           Whether to plot only the bipolar electrograms, only the
%           unipolar electrograms or both
%   'reference' 'off'|{'on'}
%           Whether to plot the reference channel, off by default userdata.electric.egm(iEgmArray,:) are selected for plotting
%           To convert from Carto point numbers to iEgmArray use
%           getIndexFromCartoPointNumber.
%
% GETEGMSATPOINTS by default returns all the electrograms of 'egmtype'. Use 
% getIndexFromCartoPointNumber to convert from point numbers to index 
% numbers.
%
% TODO: Devolve the functionality of plotOpenEPEgms to getEgmsAtPoints
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% See also PLOTELECTROGRAMS, GETINDEXFROMCARTOPOINTNUMBER,
% PLOTOPENEPEMGS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% [egmTraces, acttime, egmNames] = getEgmsAtPoint(userdata ...
%                                   , 'iegm', 1 ...
%                                   , 'egmtype', 'bip' ...
%                                   , 'reference', 'off' ...
%                                   );
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 2; % UPDATE VALUE
iEgmArray = ':';
egmtype = 'bip-uni';
reference = 'on';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch lower(varargin{i})
            case 'iegm'
                iEgmArray = varargin{i+1};
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
        case 'uni'
            for i = 1:2
                egmTraces{i,iEgm} = userdata.electric.egmUni(iEgmArray(iEgm),:,i);  
                acttime(i,iEgm) = userdata.electric.annotations.mapAnnot(iEgmArray(iEgm));
                egmNames{i,iEgm} = char(userdata.electric.names(iEgmArray(iEgm),:));
            end
        case 'bip-uni'
            i = 1;
            egmTraces{i,iEgm} = userdata.electric.egm(iEgmArray(iEgm),:);
            egmNames{i,iEgm} = char(userdata.electric.names(iEgmArray(iEgm),:));
            acttime(i,iEgm) = userdata.electric.annotations.mapAnnot(iEgmArray(iEgm));
            for i = 1:2
                egmTraces{i+1,iEgm} = userdata.electric.egmUni(iEgmArray(iEgm),:,i);  
                acttime(i+1,iEgm) = userdata.electric.annotations.mapAnnot(iEgmArray(iEgm));
                egmNames{i+1,iEgm} = char(userdata.electric.names(iEgmArray(iEgm),:));
            end
    end
    if strcmpi(reference, 'on')
        egmTraces{end+1,iEgm} = userdata.electric.egmRef(iEgmArray(iEgm),:) / 20;
        acttime(end+1,iEgm) = userdata.electric.annotations.referenceAnnot(iEgm);
        egmNames{end+1,iEgm} = 'Ref';
    end
end
egmTraces = flipud(egmTraces(:));
acttime = flipud(acttime(:));
egmNames = flipud(egmNames(:));

end