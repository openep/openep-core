function info = parse_header(header, filetype)
%parse_header Parse header data for St. Jude Export files
%   Parse the header for exported data from St. Jude export files. Written
%   to be compatible with both DxL exports and ECG exports (though the
%   latter is mostly a hack of the former)

%{
Parameters
----------
header : str
    String of the header data
filetype : 'dxl' or other
    String to identify which filetype header is being processed, as
    different exports will have different data

Returns
-------
info : struct
    Data extracted from header
%}

% Author Phil Gemmel 2020

% in some csv formats carriage return ('\r') is used as well as ('\n')
count = 0;
startIndex = 999;
while ~isempty(startIndex)
    startIndex = regexp(header,char([13 10])); % char([13 10]) = '\r\n'
    header(startIndex) = [];
    if ~isempty(startIndex)
        warning('Both carriage return and newline operators detected - carriage returns removed')
    end
    count = count + 1;
    if count > 3
        error('There is an excess number of carriage returns - please check.')
    end
end

indNL = [1, 1 + find(header==char(10))];    %#ok<*CHARTEN> %indNL = index of first character after newline
iL = 0;
startTime = NaN;
endTime = NaN;
while iL < (numel(indNL)-1)
    iL = iL + 1;
    line = header(indNL(iL):(indNL(iL+1)-2));

    % look for - "Map name:, TEXT"
    tokens = regexp(line, 'Map name\s*:\s*,\s*([a-zA-Z_.0-9\s]*)', 'tokens');
    if ~isempty(tokens)
        mapName = tokens{1}{1};
    end

    % look for - "Map type:, TEXT"
    tokens = regexp(line, 'Map type\s*:\s*,\s*(\w*)', 'tokens');
    if ~isempty(tokens)
        mapType = tokens{1}{1};
    else
        mapType = 'N/A';
    end

    % look for - "Export File Version : TEXT"
    tokens = regexp(line, 'Export File Version\s*:\s*([a-zA-Z_.0-9]*)R', 'tokens');
    if ~isempty(tokens)
        exportFileVersion = str2double(tokens{1}{1});
    end

    % look for - "St. Jude Medical. File Revision : TEXT"
    tokens = regexp(line, 'St. Jude Medical. File Revision\s*:\s*([a-zA-Z_.0-9]*)', 'tokens');
    if ~isempty(tokens)
        if exist('exportFileVersion', 'var')
            error('OPENEP/PARSE_HEADER: exportFileVersion has already been set');
        else
            exportFileVersion = str2double(tokens{1}{1});
        end
    end

    % look for - "Export Data Element : NAME"
    tokens = regexp(line, 'Export Data Element\s*:\s*(\w*)', 'tokens');
    if ~isempty(tokens)
        dataElement = tokens{1}{1};
    end

    % look for - "Export from Study : NAME"
    tokens = regexp(line, 'Export from Study\s*:\s*(\w*)', 'tokens');
    if ~isempty(tokens)
        study = tokens{1}{1};
    end

    % look for - "Export from Segment : NAME"
    %tokens = regexp(line, 'Export from Segment\s*:\s*(\w*\s*\d*)', 'tokens'); % changed 22nd May 2022
    tokens = regexp(line, 'Export from Segment\s*:\s*(.*)', 'tokens');
    if ~isempty(tokens)
        studySegment = tokens{1}{1};
    end

    % look for - "Export User Comments : TEXT"
    tokens = regexp(line, 'Export User Comments\s*:\s*(.*)', 'tokens');
    if ~isempty(tokens)
        userComments = tokens{1}{1};
    end

    % look for - "Export Start Time (h:m:s.msec) : HH:MM:SS.MS"
    tokens = regexp(line, '\s*Export Start Time \(h:m:s\.msec\)\s*:\s*(\d*:\d*:\d*\.\d*)', 'tokens');
    if ~isempty(tokens)
        startTime = tokens{1}{1};
    end

    % look for - "Export End Time (h:m:s.msec) : HH:MM:SS.MS"
    tokens = regexp(line, '\s*Export End Time \(h:m:s\.msec\)\s*:\s*(\d*:\d*:\d*\.\d*)', 'tokens');
    if ~isempty(tokens)
        endTime = tokens{1}{1};
    end

    % look for - "Export Start Time (secs usecs) : sec usec"
    tokens = regexp(line, '\s*Export Start Time \(secs usecs)\s*:\s*(\d*\s*\d*)', 'tokens');
    if ~isempty(tokens)
        if ~isempty(tokens{1}{1})
            % Extract and convert to float
            startTimeAbs = strsplit(tokens{1}{1}, ' ');
            startTimeAbs{2} = ['0.', startTimeAbs{2}];
            startTimeAbs = str2double(startTimeAbs{1})+str2double(startTimeAbs{2});
        else
            startTimeAbs = NaN;
        end
    end

    % look for - "Export End Time (secs usecs) : sec usec"
    tokens = regexp(line, '\s*Export End Time \(secs usecs)\s*:\s*(\d* \d*)', 'tokens');
    if ~isempty(tokens)
        endTimeAbs = strsplit(tokens{1}{1}, ' ');
        endTimeAbs{2} = ['0.', endTimeAbs{2}];
        endTimeAbs = str2double(endTimeAbs{1})+str2double(endTimeAbs{2});
    end

    % look for - "Total number of data points (columns): , "
    if exportFileVersion < 10.0
        [~, ind2] = regexp(line, 'Total number of data points \(columns\)\:\s*\,\s*', 'start', 'end');
        if ~isempty(ind2)
            numPoints = str2double(line(ind2+1:end));
        end
    else
        tokens = regexp(line, '# mapping pts\s*:\s*,\s*(\d*)', 'tokens');
        if ~isempty(tokens)
            numPoints = str2double(tokens{1}{1});
        end
    end
    % if numPoints still does not exist, look for number of freeze groups
    % in case this is a wave_refs file
    if ~exist('numPoints', 'var')
        tokens = regexp(line, '# freeze groups\s*:\s*,\s*(\d*)', 'tokens');
        if ~isempty(tokens)
            numPoints = str2double(tokens{1}{1});
            warning('OPENEP/PARSE_HEADER: Using number of freeze groups rather than number of points which could not be found');
        end
    end

    % look for - "This is file X of Y for map, "
    tokens = regexp(line, '\s*This is file (\d*) of (\d*) for map,(\w*)', 'tokens');
    if ~isempty(tokens)
        fileIndices = [str2double(tokens{1}{1}) str2double(tokens{1}{2})];
        mapId = tokens{1}{3};
    end

    % look for - "Seg data len:"
    if exportFileVersion < 10.0
        tokens = regexp(line, 'Seg data len\:\,(\d*.?\d*)', 'tokens');
    else
        tokens = regexp(line, 'Segment data length \(sec\)\s*:\s*,\s*(\d*.?\d*)', 'tokens');
    end
    if ~isempty(tokens)
        segmentDataLength = str2double(tokens{1}{1});
    end

    % Look for - "Exported seconds:"
    tokens = regexp(line, 'Exported seconds\:\,(\d*.?\d*)', 'tokens');
    if ~isempty(tokens)
        exportedSeconds = str2double(tokens{1}{1});
    end

    % Look for - "Sample rate:"
    tokens = regexp(line, 'Sample rate\:\,(\d*.?\d*)', 'tokens');
    if ~isempty(tokens)
        sampleFreq = str2double(tokens{1}{1});
    end

    % Look for - "CFE P-P sensitivity (mv)"
    tokens = regexp(line, 'CFE P-P sensitivity \(mv\)\,(\d*.?\d*)', 'tokens');
    if ~isempty(tokens)
        cfePpSensitivity = str2double(tokens{1}{1});
    end

    % Look for - "CFE Width (ms)"
    tokens = regexp(line, 'CFE Width \(ms\)\,(\d*.?\d*)', 'tokens');
    if ~isempty(tokens)
        cfeWidth = str2double(tokens{1}{1});
    end

    % Look for - "CFE Refractory (ms)"
    tokens = regexp(line, 'CFE Refractory \(ms\)\,(\d*.?\d*)', 'tokens');
    if ~isempty(tokens)
        cfeRefractory = str2double(tokens{1}{1});
    end

    % Look for - "Data starts in row"
    tokens = regexp(line, 'Data starts in row\s*,\s*(\d*)', 'tokens');
    if ~isempty(tokens)
        dataStartRow = str2double(tokens{1}{1});
    end
end

% if numPoints is still non-existent then throw and error
if ~exist('numPoints', 'var')
    error('OPENEP/PARSE_HEADER: Did not identify number of points or freeze groups - unknown data lenght and cannot proceed!');
end

info.exportFileVersion = exportFileVersion;
info.dataElement = dataElement;
info.study = study;
info.userComments = userComments;
info.studySegment = studySegment;
info.startTime = startTime;
info.endTime = endTime;
info.startTimeAbs = startTimeAbs;
info.endTimeAbs = endTimeAbs;
info.header = header;
if exportFileVersion >= 10.0
    info.dataStartRow = dataStartRow;
    info.mapName = mapName;
    info.mapType = mapType;
end
info.numPoints = numPoints;

% only do this next part if we have wave data, i.e. if sampleFreq is populated
if exist('sampleFreq', 'var')
    if exportFileVersion < 10.0
        info.fileIndices = fileIndices;
        info.mapId = mapId;
        info.segmentDataLength = segmentDataLength;
        info.exportedSeconds = exportedSeconds;
        info.sampleFreq = sampleFreq;
        info.CFE_PP_Sensitivity = cfePpSensitivity;
        info.CFE_Width = cfeWidth;
        info.CFE_Refractory = cfeRefractory;
    else
        info.segmentDataLength = segmentDataLength;
        info.exportedSeconds = segmentDataLength; % we no longer have exported seconds in DxL version above at least 10.0R
        info.sampleFreq = sampleFreq;
    end
end

end

