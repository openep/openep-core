function [data, labels] = read_visitag_sites(varargin)
% READ_VISITAG_SITES loads Carto3 visitag sites.txt file.
% Usage:
%   [sites] = read_visitag_sites(filename)
% Where:
%   data is a matrix of all the data from the file
%   labels is a cell array of the header lines from the file
%   filename is the full file path
%
% READ_VISITAG_SITES detailed description goes here
%
% Author: Steven Williams (2015) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------
% [data, labels] = read_visitag_sites(varargin)
% ---------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

filename = varargin{1};

if isempty(strfind(filename, 'Sites'))
    error('READ_VISITAG_SITES: filename must be a "Sites.txt" file.')
end
fid = fopen(filename, 'r');
if fid == (-1)
    error(['READ_VISITAG_SITES: Could not read the file: "' filename '"']);
end
try
    indata = textscan(fid, '%f', 'HeaderLines', 1);
    D = indata{1};
    data = reshape(D, 14, numel(D)/14)';
    fseek(fid, 0, 'bof');
    tLine = fgetl(fid);
    labels = strsplit(tLine);
    labels(strcmpi('',labels)) = [];
    fclose(fid);
catch
    fclose(fid);
    rethrow(err);
end

end