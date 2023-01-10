function info = loadprecision_shadows(filename)
%LOADPRECISION_SHADOWS Summary of this function goes here
%   Detailed explanation goes here
% If filename does not contain a file with a valid format, info is left
% empty and a warning generated with msgID = 'LoadPrecision:InvalidFile'.

    info = [];
    
    [~,~,ext] = fileparts(filename);
    if ~strcmpi(ext, '.csv')
        warning('LoadPrecision:InvalidFile','LOADPRECISION_WAVEFILE: .csv file expected');
        return
    end


    info = [];
    fileID = fopen(filename, 'r');
    if fileID == (-1); error('LOADPRECISION_SHADOWS: Could not open file.'); end
    cleanupFile = onCleanup(@()fclose(fileID));

    % read the data in manageable chunks
    maxBytes = 100000; % enough data to cover the header
    fseek(fileID, 0, 'bof');
    if maxBytes > filebytes2end(fileID); maxBytes = filebytes2end(fileID); end
    [fData, fDataSize] = fread(fileID, maxBytes, '*char'); fData = fData(1:fDataSize)';

    % do the prechecks and return if bad
    if ~local_prechecks(fData); return; end
    
    warning('Need to read the header properly. Ideally put all of this routine into (or consolidate with) loadprecision_wavefile.m?')
    info.dataElement = 'Shadows';
    info.softwareVersion = [];
    info.study = [];
    info.userComments = [];
    info.studySegment = [];
    info.startTime = [];
    info.endTime = [];
    info.header = [];
    info.sampleFreq = [];


    fseek(fileID, 0, 'bof');
    % Read each line and see if it contains something that we need
    tLine = '';
    while ~strcmpi(tLine, 'eof')
        tLine = fgetl(fileID);
        if strcmp(tLine, '')
            continue;
        else
            switch lower(tLine(1:3))
                case lower('Nam')
                    %The next lines first item should be the shadow names
                    name = {};
                    tLine = fgetl(fileID);
                    while ~strcmp(tLine, '')
                        commaPos = regexp(tLine, ',', 'once');
                        name{end+1} = regexprep(tLine(1:commaPos-1), '\s', ''); % remove spaces
                        tLine = fgetl(fileID);
                    end
                case lower('Ch1')
                    %The next lines contain the coordinates
                    tLine = fgetl(fileID);
                    iShadow = 0;
                    while ~strcmpi(tLine, 'eof')
                        iShadow = iShadow + 1;
                        numData = textscan(tLine, '%f,');
                        numData = numData{1};
                        %There are 9 numbers for each electrode:
                        %  Ch Fs Rx Ry Rz x y z Sp
                        %We are interested in:
                        %  Ch Rx Ry Rz x y z
                        for iChannel = 1:length(numData)/9
                            iBase = (iChannel-1)*9;
                            channel(iChannel, iShadow) = numData(1+iBase);
                            coordreal(iChannel, :, iShadow) = numData(3+iBase:5+iBase);
                            coordscaled(iChannel, :, iShadow) = numData(6+iBase:8+iBase);
                        end
                        tLine = fgetl(fileID);
                    end
            end
        end
    end

    for iShadow = 1:length(name)
        % Assign the shadow name
        data{iShadow}.Name = name{iShadow};

        % Assign the shadow channels
        data{iShadow}.Channels = channel(:,iShadow);

        % Find which channel entries should be removed
        tf = (coordreal(:,:,iShadow)==0);
        tf_b = (coordreal(:,:,iShadow)==0);
        if tf~=tf_b
            error('LOADVELOCITYSHADOWS: There is a discrepancy between the figures');
        end
        data{iShadow}.Channels(tf(:,1)) = [];

        cR = coordreal(:,:,iShadow);
        cR = cR(any(cR,2),:);
        data{iShadow}.RxRyRz = cR;

        cS = coordscaled(:,:,iShadow);
        cS = cS(any(cS,2),:);
        data{iShadow}.xyz = cS;
    end
    info.shadows = data;
    info.fileLoaded = filename;
end

function success = local_prechecks(fData)
    
    success = false;

    % Check Export Data Element "Export Data Element : NAME"
    tokens = regexp(fData, 'Export Data Element\s*:\s*(\w*)', 'once', 'tokens');
    goodDataElements = {'Shadows' };
    if ~isempty(tokens)
        dataElement = tokens{1};
        tf = strcmp(dataElement, goodDataElements);
        if all(not(tf))
            warning('LoadPrecision:InvalidFile','LOADPRECISION_SHADOWS: Invalid Export Data Element');
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

    if isempty(ind1) && isempty(ind2) && isempty(ind3) && isempty(ind4) && isempty(ind5)
        warning('LoadPrecision:InvalidFile','LOADPRECISION_SHADOWS: Invalid File Revision Number');
        return
    end
    
    success = true;
end
