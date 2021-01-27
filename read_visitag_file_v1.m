function [ data, header ] = read_visitag_file_v1( filepath )
% READ_VISITAG_FILE_V1 Reads a data file stored in a WiseTag or VisiTag
% directory
%
% Usage:
%   [ data, header ] = read_visitag_file_v1( filepath )
% Where:
%   filepath  - path to the file to be read
%   data  - the the file data
%   header  - the file header
%
% READ_VISITAG_FILE_V1 Reads numeric data from data files stored in WiseTag or
% VisiTag directories. These files all have a standard format - header line
% followed by data lines; tab delimited. The headers are returned in a cell
% array to allow search/indexing.
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

% read the data
data = dlmread(filepath, '', 1,0);

% read the header
fid = fopen(filepath);
line = fgetl(fid);
fclose(fid);
header = line;
header = split(header);
iRemove = false(numel(header));
for i = 1:numel(header)
    if isempty(header{i})
        iRemove(i) = true;
    end
end
header(iRemove) = [];