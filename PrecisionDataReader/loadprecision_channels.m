function [ data ] = loadprecision_channels(filename)
%LOADPRECISION_CHANNELS Summary of this function goes here
%   Detailed explanation goes here
% If filename does not contain a file with a valid format, info is left
% empty and a warning generated with msgID = 'LoadPrecision:InvalidFile'.

%   data is a m*3 matrix read from EnGuideSettings.csv. There is
%   one row per EnGuide electrode and columns [EnGuideID Channel Num] where
%   Num is the electrode number for that catheter (1 = distal).
    
    info = [];
    
    if  ~contains(filename, 'Channels' )
        warning('LoadPrecision:InvalidFile','LOADPRECISION_CHANNELS: "Channels" expected in filename.');
        return
    end

    warning('Need to read the header properly. Ideally put all of this routine into (or consolidate with) loadprecision_wavefile.m?')
    info.dataElement = 'Channels';
    info.softwareVersion = [];
    info.study = [];
    info.userComments = [];
    info.studySegment = [];
    info.startTime = [];
    info.endTime = [];
    info.header = [];
    info.sampleFreq = [];
    info.fileLoaded = filename;
    
    
    
    data = [];
    fid = fopen(filename, 'r');
    
    % Read each line and see if it contains something that we need
    tLine = '';
    while ~strcmpi(tLine, 'eof')
        tLine = fgetl(fid);
        if strcmp(tLine, '')
            continue;
        else
            switch lower(tLine(1:3))
                case lower('ID,')
                    %The next lines contain the EnGuide configurations
                    tLine = fgetl(fid);
                    ind = 0;
                    while ~strcmpi(tLine, 'eof')
                        ind = ind + 1;
                        %It is a CSV file, so scan to create a cell array
                        values = textscan(tLine, '%s', 'delimiter', ',');
                        values = values{:};
                        
                        %Remove empty cells from the end of the cell array
                        for i = length(values):-1:1
                            if ~strcmp(values(i), '')
                                break
                            else
                                values(i) = [];
                            end
                        end
                        
                        %Calculate the number of electrodes
                        numElec(ind) = (length(values) - 20) / 5;
                        
                        %Get the EnGuide number
                        id(1:numElec(ind),ind) = str2num(values{1});
                        
                        %Get the channel-number pairs
                        chstart = 21; %the column idx of the first electrode
                        chspacing = 5; %the column idx of the second electrode is chstart+chspacing
                        elecoffset = 1; %distance from chstart to electrode number
                        for i = 1:numElec(ind)
                            iCh = chstart + (i-1)*chspacing;
                            
                            ch(i,ind) = str2num(values{iCh});
                            
                            elecstr = values{iCh+elecoffset};
                            if ischar(elecstr) && strcmpi(elecstr, 'D')
                                elec = 1;
                            elseif ischar(elecstr)
                                elec = str2num(elecstr);
                            end
                            el(i,ind) = elec;
                        end
                        tLine = fgetl(fid);
                    end
            end
        end
    end
    
    id_data = id(:);
    tf = (id_data==0);
    id_data(tf) = [];
    
    ch_data = ch(:);
    ch_data(tf) = [];
    
    el_data = el(:);
    el_data(tf) = [];
    
    data = [id_data ch_data el_data];
    
    % for iShadow = 1:length(name)
    %     data{iShadow}.Name = name{iShadow};
    %
    %     data{iShadow}.Channels = channel(:,iShadow);
    %     data{iShadow}.Channels(data{iShadow}.Channels==0) = [];
    %
    %     cR = coordreal(:,:,iShadow);
    %     cR = cR(any(cR,2),:);
    %     data{iShadow}.RxRyRz = cR;
    %
    %     cS = coordscaled(:,:,iShadow);
    %     cS = cS(any(cS,2),:);
    %     data{iShadow}.xyz = cS;
    % end
    fclose(fid);
    info.channels = data;
    
end

