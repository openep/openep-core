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
% [info, varnames, data] = loadensitex_dxldata('<pathtofile>/Contact_Mapping/Map_PP_bi.csv');
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

disp(['LOADENSITEX_DXLDATA: Reading file: ' filename]);
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
[~, ind2] = regexp(fData, '****','start','end');
if isempty(ind2)
    error('End of header not found. Double check that maxBytes is large enough to cover header.')
end
indEndofHeader = ind2(end);
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
for i = 1:info.dataStartRow-1
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

% we are already at the right line in the file as we just read the header line before the data
numericColumnsToRead = tfNumHeaders;
varColumnsToRead = ~tfNumHeaders;
if isfield(info, 'mapType')
    parseMethod = 'internal'; % we are dealing with a map file
else
    parseMethod = 'regexp'; % faster for dealing with wave data
end
data = local_parsedata(fileID, varColumnsToRead, numericColumnsToRead, info.numPoints, [thisFileName ext], parseMethod);

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

function allOutput = local_parsedata(fileID, varColumnsToRead, numericColumnsToRead, nSamples, fname, parseMethod)
%   nSamples - the number of samples to read; which may be a number of
%   points or a number of freeze groups
%   columnsToRead - logical array indicating which columns will be read
%   fileID - the file ID

nNumericColToRead = sum(numericColumnsToRead);
nCol = numel(numericColumnsToRead);

maxBytes = 10 * 1024 * 1024; % read in max 10MBytes at a time
allNumericData = zeros(nSamples, nNumericColToRead, 'double');
allVarData = cell(nSamples, nCol - nNumericColToRead);   %CHANGED HERE
currentLine = 1;
remainingBytes = filebytes2end(fileID);
totalBytes = remainingBytes;
remainingData = [];
set(0,'DefaultTextInterpreter','none')
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

    switch parseMethod
        % This section needs to output reshapedData and wholeLinesRead
        % The internal method is robust and seems to work with most files
        % but is quite slow. The regexp method is much faster but fails
        % with some mapping files. 
        %
        % We need to confirm but the regexp method MIGHT work fine with all
        % wave files; in which case we will detault to using INTERNAL for
        % mapping files and REGEXP for wave files.
        case 'internal'
            % Internal method - robust but very slow
            dataChunkCellArray = parseCSVString(dataChunk);
            reshapedData = dataChunkCellArray(:, 1:nCol);  % remove extra columns

            numLinesRead = numel(reshapedData) / nCol;
            wholeLinesRead = floor(numLinesRead);

        case 'regexp'
            % split the text at commas
            dataChunkCellArray = regexp(dataChunk', ',', 'split'); % deals with successive delimiters correctly in contrast to strsplit(dataChunk', ',');

            %remove any leading or trailing empty cells if needed
            if isempty(dataChunkCellArray{1})
                dataChunkCellArray(1) = [];
            end
            if isempty(dataChunkCellArray{end})
                dataChunkCellArray(end) = [];
            end
            if strcmpi(dataChunkCellArray{end}(2:end), 'EOF')
                dataChunkCellArray(end) = [];
            end

            % % work out the valid cells
            numLinesRead = numel(dataChunkCellArray) / nCol;
            wholeLinesRead = floor(numLinesRead);

            if numLinesRead > wholeLinesRead
                % there was overhanging data, so increment nCol
                nCol = nCol + 1;
            end

            % reshape the data
            reshapedData = reshape(dataChunkCellArray(1:nCol*wholeLinesRead),[nCol, wholeLinesRead]);
            reshapedData = reshapedData';

            if numLinesRead > wholeLinesRead
                % there was overhanging data, now is the time to remove it
                reshapedData(:,end) = [];
                % and decrement nCol
                nCol = nCol-1;
            end
    end

    % Deal first with the numeric data - only keep the columns we want for signal data
    thisSignalData = reshapedData(:,numericColumnsToRead);

    % equivalent to, but much faster than, allData(currentLine:currentLine+wholeLinesRead-1,1:nColToRead) = str2double(thisEgmData);
    doubleValues = sscanf(sprintf(' %s',thisSignalData{:}),'%f',[1,Inf]);
    doubleValueReshaped = reshape(doubleValues, size(thisSignalData));
    allNumericData(currentLine:currentLine+wholeLinesRead-1,1:nNumericColToRead) = doubleValueReshaped;

    % Now deal with the variables data
    thisVarData = reshapedData(:,varColumnsToRead); %opposite of numericColumnsToRead
    allVarData(currentLine:currentLine+wholeLinesRead-1,1:(nCol) - nNumericColToRead) = thisVarData;

    % increment the current line index, waitbar and remaining bytes
    currentLine = currentLine+wholeLinesRead;
    waitbar((totalBytes-remainingBytes)/totalBytes, f);
    remainingBytes = filebytes2end(fileID);
end

% destroy the waitbar
close(f)

% assign the output
allOutput = allVarData;
if ~isempty(allNumericData) %check if we are dealing with a map or an electrogram file ...
    widthOfAllOutput = size(allOutput,2);
    for iD = 1:size(allNumericData,1)
        allOutput{iD,widthOfAllOutput+1} = allNumericData(iD,:);
    end
end

    function C = parseCSVString(s)
        % parseCSVString Parse CSV from a character vector into a cell array.
        %   C = parseCSVString(s) returns an MxN cell array of char, where each
        %   row is a CSV record and each column a field. Quoted fields and
        %   embedded commas/newlines are handled. Empty fields are preserved.
        %
        %   Input:
        %     s - character vector (single string) containing the whole CSV text.
        %
        % Example:
        %   s = 'A,"B, with comma",,C\n"D with ""quote""",E,';
        %   C = parseCSVString(s);

        if ~ischar(s) && ~isstring(s)
            error('Input must be a character vector or string.');
        end
        s = char(s);                % ensure char vector
        n = numel(s);

        rows = {};                  % cell array of rows (each row is a cell vector)
        curField = '';              % current field buffer (char)
        curRow = {};                % current row (cell array)
        inQuote = false;
        i = 1;

        while i <= n
            ch = s(i);
            if ch == '"'         % quote handling
                if inQuote
                    % possible escaped quote: lookahead
                    if i < n && s(i+1) == '"'
                        curField(end+1) = '"'; % append one quote
                        i = i + 1;             % skip the escaped quote
                    else
                        % closing quote
                        inQuote = false;
                    end
                else
                    % starting quote (enter quoted mode)
                    inQuote = true;
                end
                i = i + 1;
                continue;
            end

            if ~inQuote
                if ch == ','    % field separator
                    curRow{end+1} = curField; %#ok<AGROW>
                    curField = '';
                    i = i + 1;
                    continue;
                end

                % newline handling: support \r\n, \n, or \r
                if ch == sprintf('\r')    % CR
                    % check for CRLF
                    if i < n && s(i+1) == sprintf('\n')
                        i = i + 2;
                    else
                        i = i + 1;
                    end
                    % finish row
                    curRow{end+1} = curField; %#ok<AGROW>
                    rows{end+1,1} = curRow;    %#ok<AGROW>
                    curRow = {}; curField = '';
                    inQuote = false;
                    continue;
                elseif ch == sprintf('\n') % LF
                    i = i + 1;
                    curRow{end+1} = curField; %#ok<AGROW>
                    rows{end+1,1} = curRow;    %#ok<AGROW>
                    curRow = {}; curField = '';
                    inQuote = false;
                    continue;
                end
            end

            % normal character (either inside quotes or plain text)
            curField(end+1) = ch;
            i = i + 1;
        end

        % End of input: push remaining field/row
        % If the input ended while inside a quoted field, we treat it as finished.
        curRow{end+1} = curField;
        rows{end+1,1} = curRow;

        % Convert rows (cell of cell) into a rectangular M-by-N cell array padded with ''
        M = numel(rows);
        maxCols = 0;
        for r = 1:M
            maxCols = max(maxCols, numel(rows{r}));
        end

        C = repmat({''}, M, maxCols);
        for r = 1:M
            rowCells = rows{r};
            C(r,1:numel(rowCells)) = rowCells;
        end
    end

end