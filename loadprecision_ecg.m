function ecg = loadprecision_ecg(filename)
%loadprecision_ecg loads the ECG data
%   Loads the ECG data from the .csv file that is output by the St Jude
%   Medical software
%{
Parameters
----------
filename : str
    Name of .csv file to be read

Flags
-----

Returns
-------
ecg : struct
    ECG data, coupled with information from the original file

Revision History
----------------
File created: Philip Gemmell (2020-04-22)
%}

file_limiter = find(filename=='.');
ext = filename(file_limiter(1):end);
assert(strcmpi(ext, '.csv'), 'Expected CSV file.')

fileID = fopen(filename, 'r');
if fileID == (-1)
    error('LOADPRECISION_DXLDATA: Could not open file.')
end
cleanupFile = onCleanup(@()fclose(fileID));

% We are going to read the data in chunks from the file and then
% process it. Assume that the first chunk has enough data in it to span
% the entire 'header' region.

% read the data in manageable chunks
maxBytes = 100000; % enough data to cover the header
fseek(fileID, 0, 'bof'); % move to the beginning of the file
if maxBytes > filebytes2end(fileID)
    maxBytes = filebytes2end(fileID);
end
[fData, fDataSize] = fread(fileID, maxBytes, '*char');
fData = fData(1:fDataSize)';

% do the prechecks and return if bad
if ~loadprecision_prechecks(fData, 'ECG')
    return
end

% READ THE HEADER
% ---------------
% The 'header' finishes at the end of the last line starting with "Number of samples (rows):"
% [ind1, ~] = regexp(fData, 'Number of samples (rows)[^\n]*\n','start','end');
[ind1, ~] = regexp(fData, 'Number of samples','start','end');
if isempty(ind1)
    error('End of header not found. Double check that maxBytes is large enough to cover header.')
end
indEndofHeader = ind1(end)-2;
header = fData(1:indEndofHeader);

ecg = parse_header(header, 'ecg');

% READ THE ECG
% ------------
% ind1 = regexp(fData, 'Number of waves (columns): ,', 'start');
% ind2 = regexp(fData, 'Number of samples (rows): ,', 'start');
temp_ecg = readtable(filename);
temp_ecg(end, :) = [];
temp_ecg = removevars(temp_ecg, {'Var41'});
ecg.ecg = temp_ecg;

end

