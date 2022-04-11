function data = getAblationTimeSeriesData( userdata, iTag, varargin )
% GETABLATIONTIMESERIESDATA Is used to convert location based Visitag data
% into time series ablation data
%
% Usage:
%   data = getAblationTimeSeriesData( userdata, iTag )
% Where:
%   userdata  - an OpenEP data structure containing ablation data
%   iTag - the index of the tag of interest
%   data  - the output, by default the time indices of ablation
%
% GETABLATIONTIMESERIESDATA accepts the following parameter-value pairs
%   'type'     {'time'}|'impedance'|'ai'
%
% GETABLATIONTIMESERIESDATA is used to convert location based Visitag data
% into time series ablation data. The data in userdata.rfindex.grid is
% arranged according to the number of Visitags. So, for example, if there
% are n Visitags, stored in userdata.rfindex.tag, then there will be n
% structures in userdata.rfindex.grid. These structures contain all the
% time series data (e.g. impedance, AI) relating to an individual Visitag -
% this is the grid data.
%
% Note that the concept of grid data deals with spatial distributions of
% ablation data. So, within userdata.rfindex.grid data is primarly arranged
% according to location rather than time.
%
% For this reason, the time series data is _duplicated_ within
% userdata.rfindex.grid and needs to first be concatenated before
% duplicates are removed. We have provided the function
% 'getAblationTimeSeriesData' for this purpose.
%
% Author: Steven Williams (2022)
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

% parse input arguments
nStandardArgs = 2; % UPDATE VALUE
type = 'time';
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'type'
                type = varargin{i+1};
        end
    end
end

% store the grid data in a local variable
gridData = userdata.rfindex.grid{iTag};

% deal with the times
allTimes = [];
for i = 1:numel(gridData)
    allTimes = [allTimes; gridData(i).time]; %#ok<*AGROW>
end

% return the data; dealing with duplicate values in time
switch type
    case 'time'
        data = unique(allTimes);

    case 'impedance'
        [~, iTimes] = unique(allTimes);

        allImpedances = [];
        for i = 1:numel(gridData)
            allImpedances = [allImpedances; gridData(i).impedance];
        end
        data = allImpedances(iTimes);

    case 'ai'
        [~, iTimes] = unique(allTimes);

        allAi = [];
        for i = 1:numel(gridData)
            allAi = [allAi; gridData(i).index.value];
        end
        data = allAi(iTimes);
end

end