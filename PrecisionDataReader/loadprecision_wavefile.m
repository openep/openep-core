function info = loadprecision_wavefile(filename)
% LOADPRECISION_WAVEFILE loads a Precision wavefile (egm data / locations)
% Usage:
%   info = loadprecision_egmdata(filename)
% Where:
%   filename is the filename
%   info is the output
%
% LOADPRECISION_WAVEFILE loads Precision files that relate to catheters,
% for example electrogram data or location data for a catheter. The files
% are epected to have the following EXPORT DATA ELEMENT and otherwise an
% error will be returned. Export Data Element (not case sensitive):
%   'EP_Catheter_Bipolar_Raw', 'EPcathBIO_COMPUTED', 'EPcathBIObipol_RAW',
%   'EPcathBIObipol_FILTERED', 'EPcathBIO_RAW', 'EPcathBIO_FILTERED',
%   'ECG_RAW', 'ECG_FILTERED', 'Respiration', 'Electrode_Locations',
%   'Locations'
% If filename does not contain a file with a valid format, info is left
% empty and a warning generated with msgID = 'LoadPrecision:InvalidFile'.


% Author: Nick Linton (2009)
% Modifications
%   Steven Williams 2017    Convert to Precision
%   Nick Linton 2017        Major update
%
% Info on Code Testing:
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------


% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------
    
    info = [];
    
    [~,~,ext] = fileparts(filename);
    if ~strcmpi(ext, '.csv')
        warning('LoadPrecision:InvalidFile','LOADPRECISION_WAVEFILE: .csv file expected');
        return
    end
    
    fileID = fopen(filename, 'r');
    if fileID == (-1); error('LOADPRECISION_WAVEFILE: Could not open file.'); end
    cleanupFile = onCleanup(@()fclose(fileID));

    % We are going to read the data in chunks from the file and then process
    % it. Assume that the first chunk has enough data in it to span the entire
    % 'header' region.

    % read the data in manageable chunks
    maxBytes = 100000; % enough data to cover the header   
    fseek(fileID, 0, 'bof');
    if maxBytes > filebytes2end(fileID); maxBytes = filebytes2end(fileID); end
    [fData, fDataSize] = fread(fileID, maxBytes, '*char'); fData = fData(1:fDataSize)';
    
    % do the prechecks and return if bad
    if ~local_prechecks(fData); return; end

    % The 'header' finishes at the end of the line starting with "t_dws,t_secs,t_usecs,t_ref"
    [ind1,ind2] = regexp(fData, 't_dws,t_secs,t_usecs,t_ref,[^\n]*\n', 'once','start','end');
    if isempty(ind1)
        error('End of header not found. Double check that maxBytes is large enough to cover header.')
    end
    ind2 = ind2-1; %remove the last \n
    
    columnHeadings = fData(ind1:(ind2));
    % for the header we will include 5 lines of data to make validation
    % easier
    headerEnd = find(fData(ind2:end)==char(10),6,'first');
    header = fData(1:(ind2-1+headerEnd(end)));
    
    dataStartBytes = 1*(ind2 + 1); % 1 byte for each character.
    % fData is going to be all of the data that needs to be parsed
    fseek(fileID, dataStartBytes, 'bof');
    
    [catheters, columnsToRead, nSamples, info] = local_parseheader(header, columnHeadings);
    
    data = local_parsedata(fileID, columnsToRead, nSamples, info.dataElement);
    
    catheters = local_insertegmdata(catheters, data, columnsToRead);
    
    info.catheters = catheters;
    info.fileLoaded = filename;

end


function [catheters, columnsToRead, nSamples, info] = local_parseheader(header, columnHeadings)
    % initialise data in case it isn't present in file
    dataElement = [];
    softwareVersion = [];
    study = [];
    userComments = [];
    studySegment = [];
    startTime = [];
    endTime = [];
    
    sampleFreq = 2034.5;


    % in some csv formats carriage return ('\r' is used as well as '\n')
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
        
    indNL = 1 + find(header==char(10));    %indNL = index of first character after newline
    iL = 0;
    bipolChannels = [];
    while iL < (numel(indNL)-1)
        iL = iL + 1;
        line = header(indNL(iL):(indNL(iL+1)-2));
        
        % look for - "Export Data Element : NAME"
        tokens = regexp(line, 'Export Data Element\s*:\s*(\w*)', 'tokens');
        if ~isempty(tokens); dataElement = tokens{1}{1}; end
        
        % look for - "Export from Study : NAME"
        tokens = regexp(line, 'Export from Study\s*:\s*(\w*)', 'tokens');
        if ~isempty(tokens); study = tokens{1}{1}; end
        
        % look for - "Export from Segment : NAME"
        tokens = regexp(line, 'Export from Segment\s*:\s*(\w*)', 'tokens');
        if ~isempty(tokens); studySegment = tokens{1}{1}; end

        % look for - "Export User Comments : TEXT"
        tokens = regexp(line, 'Export User Comments\s*:\s*(\w*)', 'tokens');
        if ~isempty(tokens); userComments = tokens{1}{1}; end

        % look for - "Export from Software Version: NAME"
        tokens = regexp(line, 'Export from Software Version\s*:\s*(\w*)', 'tokens');
        if ~isempty(tokens); softwareVersion = tokens{1}{1}; end

        % look for - "Catheter[NUMBER](name, num electrodes): NAME,NUMBER"
        tokens = regexp(line, 'Catheter\[(\d*)\]\(name, num electrodes\):\s*(\w*),\s*(\d*)', 'tokens');
        if ~isempty(tokens)
            iC = str2double(tokens{1}{1}) + 1; %iC = counter for Catheter
            catheters(iC).name = tokens{1}{2};              %#ok<AGROW>
            nE = str2double(tokens{1}{3}); %nE = nElectrodes
            catheters(iC).nElectrodes = nE;                 %#ok<AGROW>
        end
        
        % look for - "used channels : \n channel,catheter name,electrode name"
        indStart = regexp(line, '\s*used channels\s*:', 'once');
        if ~isempty(indStart) %then we have the used channel list following.
            line2 = header(indNL(iL+1):(indNL(iL+2)-2));
            if strcmp('channel,catheter name,electrode name', line2)
                txtformat = '%u8 %s %s';
            elseif strcmp('channel,catheter name,electrode name,is visible', line2)
                txtformat = '%u8 %s %s %u8';
            else
                error('Unexpected text');
            end
            %read the formatted text - first find the start and end
            indStart = indNL(iL+2);
            indFinish = indStart + regexp(header(indStart:end), '\n[\s,]*\n', 'once')-1;  % find first blank line after indStart, signifies end of data
            temp = textscan(header(indStart:indFinish), txtformat, 'Delimiter', ',');
            usedChannels.channelNo = temp{1};
            usedChannels.cathName = temp{2};
            usedChannels.electrodeName = temp{3};
        end
        
        % look for - "used bipol channels : \n bipol channel,unipol channel A,unipol channel B"
        indStart = regexp(line, '\s*used bipol channels\s*:', 'once');
        if ~isempty(indStart) %then we have the used channel list following.
            line2 = header(indNL(iL+1):(indNL(iL+2)-2));
            if ~strstartcmp('bipol channel,unipol channel A,unipol channel B', line2); error('Unexpected text'); end
            %read the formatted text - first find the start and end
            indStart = indNL(iL+2);
            indFinish = indStart + regexp(header(indStart:end), '\n[\s,]*\n', 'once')-1;  % find first blank line after indStart, signifies end of data
            temp = textscan(header(indStart:indFinish), '%u8 %u8 %u8', 'Delimiter', ',');
            bipolChannels.channelNo = temp{1};
            bipolChannels.chanNoA = temp{2}; %channel number A
            bipolChannels.chanNoB = temp{3};
        end
       
        % look for - "Number of waves (columns): ,NUMBER"
        tokens = regexp(line, '\s*Number of waves \(columns\):\s*,\s*(\d*)', 'tokens');
        if ~isempty(tokens); nWaves = str2double(tokens{1}{1}); end
        
        % look for - "Number of samples (rows): ,NUMBER"
        tokens = regexp(line, '\s*Number of samples \(rows\):\s*,\s*(\d*)', 'tokens');
        if ~isempty(tokens); nSamples = str2double(tokens{1}{1}); end
        
        % look for - "Export Start Time (h:m:s.msec) : HH:MM:SS.MS"
        tokens = regexp(line, '\s*Export Start Time \(h:m:s\.msec\)\s*:\s*(\d*:\d*:\d*\.\d*)', 'tokens');
        if ~isempty(tokens); startTime = tokens{1}{1}; end

        % look for - "Export End Time (h:m:s.msec) : HH:MM:SS.MS"
        tokens = regexp(line, '\s*Export End Time \(h:m:s\.msec\)\s*:\s*(\d*:\d*:\d*\.\d*)', 'tokens');
        if ~isempty(tokens); endTime = tokens{1}{1}; end

    end
    
    % Now we need to create the column headings according to the type of
    % file that we are reading. Also adjust sampleFreq if necessary
    if strstartcmp('EPcathBIO_', dataElement)
        channels = usedChannels;
        channels.columnHeadings = cell(size(channels.channelNo));
        for i = 1:numel(channels.channelNo)
            channels.columnHeadings{i} = ['c' num2str(channels.channelNo(i))];
        end
    elseif strstartcmp('EPcathBIObipol_', dataElement) || strstartcmp('EP_Catheter_Bipolar_', dataElement)
        bipolChannels.columnHeadings = cell(size(bipolChannels.channelNo));
        bipolChannels.cathName = cell(size(bipolChannels.channelNo));
        for i = 1:numel(bipolChannels.channelNo)
            a = find(bipolChannels.chanNoA(i) == usedChannels.channelNo, 1, 'first');
            b = find(bipolChannels.chanNoB(i) == usedChannels.channelNo, 1, 'first');
            
            nameA = usedChannels.cathName{a};
            nameB = usedChannels.cathName{b};
            elecA = usedChannels.electrodeName{a};
            elecB = usedChannels.electrodeName{b};
            
            if ~strcmp(nameA,nameB)
                error('Code not written to read bipoles across different catheters.')
            end
            bipolChannels.cathName{i} = nameA;
            bipolChannels.columnHeadings{i} = [nameA '(' elecA '-' elecB ')_c' num2str(bipolChannels.channelNo(i))];
        end
        channels = bipolChannels;
    elseif strstartcmp('ECG_', dataElement)
        % create a 'catheter' for the 12 lead ECG
        catheters.name = 'ECG';
        catheters.nElectrodes = 12;
        channels.cathName = repmat({'ECG'},12,1);
        channels.electrodeName = {'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'}';
        channels.columnHeadings = {'I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6'}';
    elseif strstartcmp('Locations', dataElement) || strstartcmp('Electrode_Locations', dataElement)
        channels = usedChannels;
        channels.columnHeadings = cell(numel(channels.channelNo),3);
        for i = 1:numel(channels.channelNo)
            channels.columnHeadings{i,1} = ['c' num2str(channels.channelNo(i)) 'x'];
            channels.columnHeadings{i,2} = ['c' num2str(channels.channelNo(i)) 'y'];
            channels.columnHeadings{i,3} = ['c' num2str(channels.channelNo(i)) 'z'];
        end
        sampleFreq = 2034.5/20;
    elseif strstartcmp('Respiration', dataElement)
        catheters.name = 'Respiration';
        catheters.nElectrodes = 2;
        channels.cathName = repmat({'Respiration'},2,1);
        channels.electrodeName = {'absolute','percentage'}';
        channels.columnHeadings = {'absolute','percentage'}';
        sampleFreq = 2034.5/20;
    else
        error(['"' dataElement '"' ' not found.'])
    end
    
    % Next we will work out the column number for each of the column
    % headings that we want.
    commaPositions = (columnHeadings == ',');
    columnNumbers = 1 + cumsum(commaPositions);
    columnsToRead = false(1,nWaves);
    
    for iC = 1:numel(catheters)
        name = catheters(iC).name;
        indChannels = strcmp(name, channels.cathName);
        
        catheters(iC).egmNames = channels.columnHeadings(indChannels,:);      
        catheters(iC).egmColumnNumbers = NaN(size(catheters(iC).egmNames));                     
        
        for iE = 1:numel(catheters(iC).egmColumnNumbers)
            txt = [',' catheters(iC).egmNames{iE} ','];
            indStart = regexp(columnHeadings, regexptranslate('wildcard', txt), 'once');
            
            if isempty(indStart)
                beep
                warning(['Column called' txt ' was not found']);
            else
                catheters(iC).egmColumnNumbers(iE) = columnNumbers(indStart);
                columnsToRead(columnNumbers(indStart)) = true;
            end
        end
    end

    info.dataElement = dataElement;
    info.softwareVersion = softwareVersion;
    info.study = study;
    info.userComments = userComments;
    info.studySegment = studySegment;
    info.startTime = startTime;
    info.endTime = endTime;
    info.header = header;
    info.sampleFreq = sampleFreq;
    
end
   
function allData = local_parsedata(fileID, columnsToRead, nSamples, dataElement)
    % The first column contains the time in hh:mm:ss.ms format, all other
    % columns can be read as %f format. In order to use sscanf, we need to
    % generate the correct FORMAT for each row.
    if columnsToRead(1); error('Code not written to read first column of data.'); end
    formatFirstCol = uint8('%*u:%*u:%*f,'); %this reads and discards the first column
    
    columnRead = uint8('%f,');
    columnLeave = uint8('%*f,');
    
    nColToRead = sum(columnsToRead);
    nColToLeave = numel(columnsToRead) - nColToRead - 1; % -1 because we are dealing with first column separately
    
    totalLength = numel(formatFirstCol) + nColToRead*numel(columnRead) + nColToLeave*numel(columnLeave);
    format = zeros(1,totalLength, 'uint8');
    
    format(1:numel(formatFirstCol)) = formatFirstCol;
    lastPlace = numel(formatFirstCol);
    
    for i = 2:numel(columnsToRead)
        if columnsToRead(i)
            format(lastPlace + (1:numel(columnRead))) = columnRead;
            lastPlace = lastPlace + numel(columnRead);
        else %don't read
            format(lastPlace + (1:numel(columnLeave))) = columnLeave;
            lastPlace = lastPlace + numel(columnLeave);
        end
    end
    
    format = char(format);
    
    % read in data in chunks
    allBytesToRead = filebytes2end(fileID);
    % read in max 10MBytes at a time
    maxBytes = 10 * 1024 * 1024;
    
    remainder = [];
    
    allData = zeros(nColToRead*nSamples,1,'single');
    lastDataPosition = 0;
    
    hBar = waitbar(0,['Reading data element: "' dataElement '"']);
    set(findall(hBar, 'type', 'text'), 'Interpreter', 'none')
    cleanupWaitBar = onCleanup(@()close(hBar));
    
    while~feof(fileID)
        remainingBytes = filebytes2end(fileID);
        waitbar(1-remainingBytes/allBytesToRead, hBar);
        
        bytesToRead = min([maxBytes, remainingBytes+1]); %the +1 ensures we read into the end of the file.
        fileText = [  remainder ; fread(fileID, bytesToRead, '*char')  ];
        
        endOfLastLine = numel(fileText);
        while fileText(endOfLastLine) ~= char(10)
            endOfLastLine = endOfLastLine - 1;
        end
        
        remainder = fileText((endOfLastLine+1):end); %work out hanging line
        fileText = fileText(1:(endOfLastLine-1)); %truncate to last full line
        
        [data, count] = sscanf(fileText, format);
        
        if rem(count,nColToRead) ~=0
            error('Unexpected amount of data read.')
        end
        
        allData(lastDataPosition + (1:count)) = data;
        lastDataPosition = lastDataPosition + count;
        
        clear data
    end
    
    if lastDataPosition ~= nColToRead*nSamples
        beep()
        warning('Not all samples read.')
    end
    
    % fread and sscanf reads the data into a long column (not row)
    
    allData = reshape(allData, nColToRead, nSamples);
    allData = allData';

end
    
function catheters = local_insertegmdata(catheters, data, columnsToRead)
% data has only got egms for columnsToRead == true
    nSamples = size(data,1);
    old2new = cumsum(columnsToRead);
    for iC = 1:numel(catheters)
        newColNumbers = old2new(catheters(iC).egmColumnNumbers);
        if ~isvector(newColNumbers)
            egms = double(data(:,newColNumbers));
            egms = reshape(egms, [nSamples , size(newColNumbers)]);
            catheters(iC).egms = single(egms);
        else
            catheters(iC).egms = data(:,newColNumbers);
        end
    end
end

function success = local_prechecks(fData)
    
    success = false;
    
    % Check Export Data Element "Export Data Element : NAME"
    tokens = regexp(fData, 'Export Data Element\s*:\s*(\w*)', 'once', 'tokens');
    goodDataElements = {'EP_Catheter_Bipolar_Raw', 'EPcathBIO_COMPUTED', 'EPcathBIObipol_RAW', 'EPcathBIObipol_FILTERED', 'EPcathBIO_RAW', 'EPcathBIO_FILTERED', 'ECG_RAW', 'ECG_FILTERED', 'Respiration', 'Electrode_Locations', 'Locations' };
    if ~isempty(tokens)
        dataElement = tokens{1};
        tf = strcmp(dataElement, goodDataElements);
        if all(not(tf))
            warning('LoadPrecision:InvalidFile','LOADPRECISION_WAVEFILE: Invalid Export Data Element');
            return
        end
    end
    
    % Check file version
    indNewLine = find(fData==char(10),1); %#ok<CHARTEN>
    firstLine = fData(1:(indNewLine-1));
    ind1 = regexp(firstLine, 'St\.\s*Jude Medical\.\s*File Revision\s*:\s*4\.1\.research', 'once');
    ind2 = regexp(firstLine, 'St\.\s*Jude Medical\.\s*File Revision\s*:\s*5\.0\.1', 'once');
    ind3 = regexp(firstLine, 'St\.\s*Jude Medical\.\s*File Revision\s*:\s*5\.2', 'once');
    ind4 = regexp(firstLine, 'Export\s*File\s*Version\s*:\s*5\.0R', 'once');
    ind5 = regexp(firstLine, 'Export\s*File\s*Version\s*:\s*5\.1R', 'once');
    ind6 = regexp(firstLine, 'St\.\s*Jude Medical\.\s*File Revision\s*:\s*5\.6', 'once');

    
    if isempty(ind1) && isempty(ind2) && isempty(ind3) && isempty(ind4) && isempty(ind5) && isempty(ind6)
        warning('LoadPrecision:InvalidFile','LOADPRECISION_WAVEFILE: Invalid File Revision Number');
        return
    end
    
    success = true;
    
end