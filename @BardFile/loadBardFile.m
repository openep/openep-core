function loadBardFile(hB)
    
    if ~isempty(hB.ChDataFileMap)
        error('@BARDFILE/LOADBARDFILE: there is already a datafilemap - write code to delete it if you want to do this.')
    end
    
    [info, channel, dataPosition] = getbarddetails(hB);

    hB.StartTime = info.tStart;
    hB.NChannels = info.nChannels;
    hB.NSamples = info.nSamples;
    hB.SampleRate = info.sampleRate;

    hB.ChRange = zeros(1,hB.NChannels);
    hB.ChLow = zeros(1,hB.NChannels);
    hB.ChHigh = zeros(1,hB.NChannels);

    for i = 1:hB.NChannels
        hB.ChName{i} = channel(i).label;
        hB.ChRange(i) = channel(i).range;
        hB.ChLow(i) = channel(i).low;
        hB.ChHigh(i) = channel(i).high;
        if channel(i).sampleRate ~= hB.SampleRate
            error('@BARDFILE/BARDFILE: Not all channels have the same sample rate!')
        end
    end
    
    createEmptyMemMap(hB);
    
    %read the data into the memory map file in manageable chunks
    if strcmp(computer, 'PCWIN') || strcmp(computer, 'PCWIN64')
        [userview, ~] = memory;
        maxBytes = floor( userview.MaxPossibleArrayBytes / 10 );   %10 to allow for safety factor
    elseif strcmp(computer, 'GLNX86') || strcmp(computer, 'GLNXA64')
        [~, w] = unix('vmstat');
        stats = regexp(w, '[0-9]*', 'match');
        maxBytes = floor( (str2double(stats{4})*4096) / 10 );
    elseif strcmp(computer, 'MACI64')
        [~, w] = unix('vm_stat | grep Pages');
        stats = regexp(w, '[0-9]*', 'match');
        maxBytes = floor( (str2double(stats{1})*4096) / 10 );
    else
        error('LOADBARDFILE: Unrecognised platform');
    end
    
    s = dir(hB.FileName);
    filesize = s.bytes;
    
    if maxBytes>filesize
        maxBytes = filesize;
    end

    
    fid = fopen(hB.FileName);
    fseek(fid, dataPosition, 'bof');
    
    fDataRemainder = [];
    linesEntered = 0;
        
    while ~feof(fid)
        
        [fData, fDataSize] = fread(fid, maxBytes, '*char');
        fData = fData(1:fDataSize)';
        fData = [fDataRemainder fData]; %#ok<AGROW>
        
        newline = char(10);
        endLastLine = find(fData == newline, 1, 'last');
        
        fDataRemainder = fData(endLastLine:end);
        
        fData(fData == ',') = ' ';  % for sscanf, replace commas with spaces
        data = sscanf( fData(1:endLastLine), '%ld');   %read it into int64 (%ld) - easier to convert to int16
        data = int16(data);
        
        nLines = numel(data) / hB.NChannels;
        if nLines-floor(nLines) > 0
            error('@BARDFILE\LOADBARDFILE: wrong right dimensions of data.')
        end
        
        data = reshape(data, hB.NChannels, nLines)';
        hB.ChDataFileMap.Data.a2d( linesEntered + (1:nLines), : ) = data;
        linesEntered = linesEntered + nLines;
    end

    fclose(fid);
    clear data %data is a large array so get rid of it
    
    %%%%%%%%%%%%%%%%%%%%
    if strcmpi(info.chPaced, 'no entry in file') && isempty(hB.ChStim)
        chList = {'none' hB.ChName{:}}; %#ok<CCAT>
        [~, fNameShort, ~] = fileparts(hB.FileName);
        [selection,ok] = listdlg(     'ListString' , chList ...
            , 'SelectionMode' , 'single' ...
            , 'ListSize' , [300, 600] ...
            , 'InitialValue' , 1 ...
            , 'PromptString', ['No pacing channel in: - ' fNameShort]  ...
            , 'OKString' , 'Add selection to txt file' ...
            , 'CancelString' , 'Leave original file alone' ...
            );
        if ok
            if selection == 1
                hB.PrivateChStim = NaN;
            else
                hB.PrivateChStim = selection - 1;
            end
            txt = chList{selection};
            txt = ['Channel paced: ' txt];
            inserttextintotextfile(hB.FileName, 4, txt);
        end
    elseif strcmpi(info.chPaced, 'no entry in file')
        warning('The BardFile object has a chPaced, but there isn''t one in the Bard text file. chPaced will be used.')
    else
        hB.ChStim = info.chPaced;
    end
    
    %%%%%%%%%%%%%%%%%%%%



end


