function [info, channel, dataPosition] = getbarddetails(obj)
% LOADBARDTEXT loads an exported Bard text file into variables.
% Usage:
%   [header, channelDetails, data] = loadbardtext(filename)
% Where:
%   info is a structure
%       .chPaced is the paced channel
%       .sampleRate is the sample rate in Hz
%       .nChannels is the number of channels
%       .nSamples is the number of samples in each channel
%       .tStart is the start time (in text format)
%   channel is an array of structures, length - info.nChannels
%       .sampleRate is the sample rate in Hz
%       .label is the title given to the channel
%   dataPosition is the file position of the A/D data
%
% LOADBARDTEXT loads the basic data - it could be improved in the future in
% order to obtain further info.

% Author: Nick Linton (2009)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------


%use warnings with the following format
%disp('MYFUNCTION: warning.')

filename = obj.FileName;

fid = fopen(filename);
if fid == (-1)
    error('LOADBARDTEXT: the file was not found.')
end
    
tLine = '';
while not(strcmp('[Header]', tLine))
    tLine = fgetl(fid);
    if tLine == -1
        error('LOADBARDTEXT: the file header was not found.')
    end
end
    


% Read the header.
% The Header is preceded by the text [Header] on a separate line:
%     [Header]
%     File Type: 1
%     Version: 1
%     Channel paced: HRAd   %% this is an optional extra that can be added manually
%     Channels exported: 6
%     Samples per channel: 41824
%     Start time: 16:31:47
%     End time: 16:32:29
%     Ch. Info. Pointer: 320
%     Stamp Data: T
%     Mux format: 0
%     Mux Block Size: 
%     Data Format 1
%     Sample Rate: 1000Hz
%     Channel #:   1                %start of channel 1 header
%     Label: I
%     Range: 5mv 
%     Low: .05Hz
%     High: 25Hz
%     Sample rate: 1000Hz
%     Color: 00FFFF
%     Scale: -6
%     Channel #:   2                %etc

% Main Details block
mainDetailsText =   { ...
                        'File Type:' ...
                        'Version:' ...
                        'Channel paced:' ...
                        'Channels exported:' ...
                        'Samples per channel:' ...
                        'Start time:' ...
                        'End time:' ...
                        'Ch. Info. Pointer:' ...
                        'Stamp Data:' ...
                        'Mux format:' ...
                        'Mux Block Size:' ...
                        'Data Format 1' ...
                        'Sample Rate:' ...
                        'Channel #:' ...
                    };

reachedChannel = false;

info = struct('sampleRate', 0, 'nChannels',0 , 'nSamples',0 , 'tStart',0 , 'chPaced','');
info.chPaced = 'no entry in file';
while reachedChannel == false
    textLine = fgetl(fid);
    for i = 1:length(mainDetailsText)
        isMatch = strstartcmpi(mainDetailsText{i}, textLine);
        if isMatch
            switch lower(mainDetailsText{i})
                case 'channel paced:'
                    info.chPaced = strtrim(textLine(15:end));
                case 'sample rate:'
                    info.sampleRate = str2double(strtrim(textLine(14:(end-2))));
                case 'channels exported:'
                    info.nChannels = str2double(strtrim(textLine(19:end)));
                case 'samples per channel:'
                    info.nSamples = str2double(strtrim(textLine(21:end)));
                case 'start time:'
                    info.tStart = (strtrim(textLine(12:end)));
                case 'channel #:'
                    reachedChannel = true;   %This will exit the while loop.
                % I haven't added in all of the info - not needed at the moment
                % but may be useful in the future.
            end
        end
    end
end

% Now the channel info
channelDetailsText = { ...
                        'Channel #:' ...
                        'Label:' ...
                        'Range:' ... 
                        'Low:' ...
                        'High:' ...
                        'Sample rate:' ...
                        'Color:' ...
                        'Scale:' ...
                        '[Data]' ...
                    };

channelNumber = 1;
reachedData = false;

channel(1,info.nChannels) = struct('label','' , 'range',0 , 'low',0 , 'high',0 , 'sampleRate',0 );

while reachedData == false
    textLine = fgetl(fid);
    for i = 1:length(channelDetailsText)
        isMatch = strstartcmpi(channelDetailsText{i}, textLine);
        if isMatch
            switch lower(channelDetailsText{i})
                case 'label:'
                    channel(channelNumber).label = strtrim(textLine(7:end));
                case 'range:'
                    channel(channelNumber).range = str2double(textLine(7:(end-3)));
                case 'low:'
                    channel(channelNumber).low = str2double(textLine(5:(end-2)));
                case 'high:'
                    if textLine(end-2) == 'k' %then the units are kHz
                        channel(channelNumber).high = 1000* str2double(textLine(6:(end-3)));
                    else
                        channel(channelNumber).high = str2double(textLine(6:(end-2)));
                    end
                case 'sample rate:'
                    channel(channelNumber).sampleRate = str2double(textLine(13:(end-2)));
                case 'channel #:'
                    channelNumber = channelNumber + 1;
                case '[data]'
                    reachedData = true;     %This will exit the while loop.
                % I haven't added in all of the info - not needed at the moment
                % but may be useful in the future.
            end
        end
    end
end
dataPosition = ftell(fid);
fclose(fid);
            
            
            