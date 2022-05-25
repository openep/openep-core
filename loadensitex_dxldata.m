function [info, varnames, data] = loadensitex_dxldata(filename)
% LOADPRECISION_DXLDATA loads the map stored in an EnSiteX DxL file.
% Usage:
%   [info, points, egms] = loadprecision_dxldata(filename)
% Where:
%   filename is the filename
%   info - contains header information about the file
%   varnames - contains the variable names stored in the data
%   data - contains the data
%
% LOADENSITE_DXLDATA Columns in egms correspond to indices in points
%
% Author: Steven Williams (2022)
% Modifications -

% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

info = [];
varnames = [];
data = [];

[~, thisFileName, ext] = fileparts(filename);
if ~strcmpi(ext, '.csv')
    warning('LoadPrecision:InvalidFile',...
        'LOADENSITEX_DXLDATA: .csv file expected');
    return
end

fileID = fopen(filename, 'r');
if fileID == (-1)
    error('LOADENSITEX_DXLDATA: Could not open file.')
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
if ~loadensitex_prechecks(fData, 'DxL')
    return
end

% READ THE HEADER
% ---------------
% The 'header' finishes at the end of the last line starting with "****,"
[ind1, ~] = regexp(fData, '****','start','end');
if isempty(ind1)
    error('End of header not found. Double check that maxBytes is large enough to cover header.')
end
indEndofHeader = ind1(end);
header = fData(1:indEndofHeader);

% Parse the header
info = parse_header(header, 'dxl');
info.filename = filename;

% % This section commented out since map data no longer appears to exist
% % READ THE MAP DATA
% % -----------------
% % The map data starts at the line that starts 'pt number:' and ends before
% % the line that starts 'Seg data len:';
% ind1 = regexp(fData, 'pt number:', 'start');
% ind2 = regexp(fData, 'Seg data len:', 'start');
%
% if isempty(ind1) || isempty(ind2)
%     points = [];
% else
%     % Parse the map data: this section is left over from Precision export
%     % formats and may not work for EnSiteX export formats
%     mapdatastring = fData(ind1:ind2-2);
%     points = local_parsemapdata(mapdatastring, info.numPoints);
% end


% READ THE ELECTROGRAMS / OTHER DATA
% ----------------------------------

% Read the header line at info.dataStartRow
fseek(fileID, 0, 'bof');
for i = 1:info.dataStartRow-2
    fgetl(fileID);
end
dataHeaderRowLine = fgetl(fileID);
dataHeaders = regexp(dataHeaderRowLine,',','split'); % previously, this was: dataHeaders = strsplit(dataHeaderRowLine, ',');

% Tidy up the heading data
if strcmpi(dataHeaderRowLine(end), ',')
    dataHeaders(end) = [];
end
if strcmpi(dataHeaders(end), '...')
    dataHeaders(end) = [];
end

% If wave data exists; i.e. if sample frequency has been set, we will read
% the wave data separately
if isfield(info, 'sampleFreq')
    headersAsNumbers = str2double(dataHeaders);
    tfNumHeaders = ~isnan(headersAsNumbers);
else
    tfNumHeaders = false(size(dataHeaders));
end

% Return the variable headings, concatenating the ending numeric data into
% a single variable. Note that we currently are assuming that NONE of the
% other header names will be numeric, this may not always be the case.
if isfield(info, 'sampleFreq')
    varnames = dataHeaders(~tfNumHeaders);
    varnames{end+1} = 'signals';
else
    varnames = dataHeaders;
end

numericColumnsToRead = tfNumHeaders;
varColumnsToRead = ~tfNumHeaders;
% we are already at the right line in the file as we just read the header line before the data

data = local_parsedata(fileID, varColumnsToRead, numericColumnsToRead, info.numPts, [thisFileName ext]);

end


function points = local_parsemapdata(mapdatastring, numPoints)
% create cell data arrays of fieldnames and data
C = textscan(mapdatastring, repmat('%s',1,numPoints+1), 'delimiter', ',', 'CollectOutput', true);
C = C{1};
fieldNames = C(:,1);
fieldNames = regexprep(fieldNames,'[^\w'']',''); %remove whitespace and punctuation
rawdata = C(:,2:end);
rawdata([1, 7:23, 25, 27:28],:) = cellfun(@(s) {str2double(s)},rawdata([1, 7:23, 25, 27:28],:)); %convert strings to doubles

% create the points structure
for iPt=1:size(rawdata,2)
    for iFn=1:length(fieldNames)
        points(iPt).(fieldNames{iFn}) = rawdata{iFn,iPt};
    end
end
end

function allOutput = local_parsedata(fileID, varColumnsToRead, numericColumnsToRead, nSamples, fname)
%   nSamples - the number of samples to read; which may be a number of
%   points or a number of freeze groups
%   columnsToRead - logical array indicating which columns will be read
%   fileID - the file ID

nNumericColToRead = sum(numericColumnsToRead);
nCol = numel(numericColumnsToRead);

maxBytes = 10 * 1024 * 1024; % read in max 10MBytes at a time
allNumericData = zeros(nSamples, nNumericColToRead, 'double');
allVarData = cell(nSamples, nCol - nNumericColToRead);
currentLine = 1;
remainingBytes = filebytes2end(fileID);
totalBytes = remainingBytes;
remainingData = [];
f = waitbar(0, ['Loading data from file: ' fname]);
while remainingBytes>0
    % Read chunk of data
    bytesToRead = min([maxBytes, remainingBytes+1]);     % The +1 ensures we read into the end of the file.
    dataChunk = fread(fileID, bytesToRead, '*char');

    % Find the final newline character
    iNewLine = regexp(dataChunk', '\n');
    lastNewLine = iNewLine(end);

    % identify the remaining data for the next time round
    temp = dataChunk(lastNewLine+1:end);

    % removing overhanging data (since the data chunk will not be an exact number of lines)
    dataChunk(lastNewLine:end) = [];

    % add on remaining data from last time if appropriate
    if ~isempty(remainingData)
        dataChunk = [remainingData; dataChunk]; %#ok<*AGROW>
    end

    % save the remaining data for the next time round
    remainingData = temp;

    % split the text at commas
    dataChunkCellArray = regexp(dataChunk', ',', 'split'); % deals with successive delimiters correctly in contrast to strsplit(dataChunk', ',');

    % remove any leading or trailing empty cells if needed
    if isempty(dataChunkCellArray{1})
        dataChunkCellArray(1) = [];
    end
    if isempty(dataChunkCellArray{end})
        dataChunkCellArray(end) = [];
    end
    if strcmpi(dataChunkCellArray{end}(2:end), 'EOF')
        dataChunkCellArray(end) = [];
    end

    % work out the valid cells
    numLinesRead = numel(dataChunkCellArray) / nCol;
    wholeLinesRead = floor(numLinesRead);

    % check if we need to insert extra cells

    % reshape the data
    reshapedData = reshape(dataChunkCellArray,[nCol, wholeLinesRead]);
    reshapedData = reshapedData';

    % Deal first with the numeric data -----

    % only keep the columns we want for signal data
    thisSignalData = reshapedData(:,numericColumnsToRead);

    % equivalent to, but much faster than
    %allData(currentLine:currentLine+wholeLinesRead-1,1:nColToRead) = str2double(thisEgmData);
    doubleValues = sscanf(sprintf(' %s',thisSignalData{:}),'%f',[1,Inf]);
    doubleValueReshaped = reshape(doubleValues, size(thisSignalData));
    allNumericData(currentLine:currentLine+wholeLinesRead-1,1:nNumericColToRead) = doubleValueReshaped;

    % Now deal with the variables data -----

    thisVarData = reshapedData(:,~numericColumnsToRead);
    allVarData(currentLine:currentLine+wholeLinesRead-1,1:nCol - nNumericColToRead) = thisVarData;

    % increment the current line index, waitbar and remaining bytes
    currentLine = currentLine+wholeLinesRead;
    waitbar((totalBytes-remainingBytes)/totalBytes, f);
    remainingBytes = filebytes2end(fileID);
end

% destroy the waitbar
close(f)

% assign the output
allOutput = allVarData;
widthOfAllOutput = size(allOutput,2);
for i = 1:size(allNumericData,1)
    allOutput{i,widthOfAllOutput+1} = allNumericData(i,:);
end

end