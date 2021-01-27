function [ newData ] = scaleData( data, limits, varargin )
%     **********************************************************************
%     * The contents of this package are copyright (c) King's College London
%     *
%     * They may not be reproduced, distributed, modified or sold for any 
%     * purpose
%     *
%     * Author: 	Dr S. E. Williams
%     * Address: 	Division of Imaging Sciences and Biomedical Engineering,
%     *          	King's College London
%     *			St Thomas' Hospital
%     *
%     * March 2014
%     **********************************************************************
%SCALEDATA performs linear scaling of a data series to the desired limits.
% Usage:
%   newData = scaleData(data, limits)
% Where:
%   data - is the input data series
%   limits - the desired bounds of the rescaled data series, [lower upper]
%
% SCALEDATA Accepts the following parameter-value pairs:
%   'scalingperiod' -   [a:b] indexes into data(a:b) to provide the range used for
%                       scaling the data. If e.g. the data has artefacts in
%                       it, then setting scalingperiod to something
%                       appropriate can avoid those artefacts
%   'clip'  - {true}|false if a custom scalingperiod has been specified
%             then setting 'clip' to true will clip any data not within
%             limits
%
% Author: Steven Williams (2012)
% Modifications -
%
% Info on Code Testing:
% ---------------------
% test code
% ---------------------
%
% Use warnings with the following format
%    disp('SCALEDATA: warning.')
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% parse the input arguments
sp = [];
clip = true;
if nargin>2
    for i = 1:2:numel(varargin)
       switch varargin{i}
           case 'scalingperiod'
               sp = varargin{i+1};
           case 'clip'
               clip = varargin{i+1};
       end
    end
end

% if the data has no range, return the data at max or min of limits
if min(data) == max(data)
    newData = data;
    if nanmin(data)>=limits(2)
        newData(~isnan(newData)) = limits(2);
    elseif nanmax(data)<=limits(1)
        newData(~isnan(newData)) = limits(1); 
    end
    return;
end

% Work out max and min
if isempty(sp)
    mx = max(data(:));
    mn = min(data(:));
else
    dataTemp = data(:);
    mx = max(dataTemp(sp));
    mn = min(dataTemp(sp));
end

% Scale the data onto [0 1]
n = (data - mn) / (mx - mn);

% Scale the data onto range(limits)
n = n * range(limits);

% Start the data at limits(1)
newData = n + limits(1);

% Clip the data
if clip
   newData(newData>limits(2)) = limits(2);
   newData(newData<limits(1)) = limits(1);
end

end

