function varargout = read_ecgfile_v4(varargin)
% READ_ECGFILE loads this Carto3 ecg file.
% Usage:
%   [headerInfo channelVoltages] = read_ecgfile(filename)
% Where:
%   headerInfo - a structure with the following fields ...
%       .uniMapChannel      - Unipolar Mapping Channel
%       .uniMapChannel2     - Unipolar Mapping Channel 2
%       .bipMapChannel      - Bipolar Mapping Channel
%       .refChannel         - Reference Channel
%       .channelNames       - names of each voltage channel e.g. 'CS1-CS2'
%       .channelNamesFull   - include the pin number e.g. 'CS1-CS2(101)'                 
%       .gain               - voltage_data * gain = voltage
%       .nSamples           - 2500 in current versions
%  channelVoltages - 2500*n integer array, where n is number of channels.
%                    See comments on channelNames and gain above.

% Author: Nick Linton (2013) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 
% 2023 - NL - changed to provide all of the raw information so that it can
% be stored in integer format (i.e. this file no longer multiplies by the
% gain).

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------


% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

headerInfo.nSamples = 2500;

filename = varargin{1};
fid = fopen(filename, 'r', 'ieee-le', 'UTF-8'); % MUCH faster than fid = fopen(filename, 'r')
if fid == (-1)
    error(['READ_ECGFILE: Could not read the file: "' filename '"']);
end
try
    line1 = fgetl(fid);

    % line 1
    line1(isspace(line1)) = [];
    if ~startsWith(line1,'ECG_Export_4.0','IgnoreCase',true) && ~startsWith(line1,'ECG_Export_4.1','IgnoreCase',true)
        error('READ_ECGFILE: The version number in the txt file is unexpected.') %#ok<*WNTAG>
    end

    line2 = fgetl(fid);
    line3 = fgetl(fid);

    % if line 
    if startsWith(line1,'ECG_Export_4.0','IgnoreCase',true)
        line4 = fgetl(fid);
    end
    if nargout>=2
        vData = fread(fid,'*char')';
    end
    fclose(fid);
catch err
    fclose(fid);
    rethrow(err)
end
    
% Check lines have expected information and retrieve it (we have already
% checked line 1)

% line 2
line2(isspace(line2)) = [];
if startsWith(line2,'rawecgtomv(gain)=0.003000','IgnoreCase',true)
    headerInfo.gain = 0.003;
else
    error('READ_ECGFILE: Unexpected statement about gain.') %#ok<*WNTAG>
end

% line 3 (optional)
if startsWith(line1,'ECG_Export_4.0','IgnoreCase',true)
    pattern = [ 'Unipolar Mapping Channel=(?<uni>\S*)\s*', ...
        'Bipolar Mapping Channel=(?<bip>\S*)\s*', ...
        'Reference Channel=(?<ref>\S*)\s*'];
    result = regexpi(line3,pattern,'names');
    if ~isempty(result)
        headerInfo.uniMapChannel = result.uni;
        headerInfo.bipMapChannel = result.bip;
        headerInfo.refChannel = result.ref;

        headerInfo.uniMapChannel2 = incrementUnipoleName(headerInfo.uniMapChannel);

        % choose the correct next line
        nextLine = line4;
    end
elseif startsWith(line1,'ECG_Export_4.1','IgnoreCase',true)
    nextLine = line3;
end

% line 4 (or line 3)
[headerInfo.channelNamesFull, nomatch] = regexpi(nextLine,'([\w-]*\(\d*\))','match','split');
% check that nomatch strings are only white space characters
test = regexp(nomatch,'\S*');
for i = 1:numel(test)
    if ~isempty(test{i}); error('Some text in channel names was missed.'); end
end
headerInfo.channelNamesFull = headerInfo.channelNamesFull';
% remove the connection number (if that's what it is),
%                                          e.g. 'CS1-CS2(101)' -> 'CS1-CS2'
headerInfo.channelNames = cell(size(headerInfo.channelNamesFull));
for i = 1:numel(headerInfo.channelNamesFull)
    result = regexp(headerInfo.channelNamesFull{i}, '(?<shortName>.*)\(\d*\)','names');
    headerInfo.channelNames{i} = result.shortName;
end

varargout{1} = headerInfo;
if nargout<2
    return
end

% the rest is the voltages
nChannels = numel(headerInfo.channelNames);
nSamples = headerInfo.nSamples;
voltages = sscanf(vData, '%d');
if numel(voltages) ~= nChannels * nSamples
    error('READ_ECGFILE: Unexpected size of voltage data.')
end
voltages = reshape(voltages, nChannels, nSamples)';

varargout{2} = voltages;

