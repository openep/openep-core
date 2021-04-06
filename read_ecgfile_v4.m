function varargout = read_ecgfile_v4(varargin)
% READ_ECGFILE loads this Carto3 ecg file.
% Usage:
%   channelNames = read_ecgfile(filename);
%   [channelNames channelVoltages] = read_ecgfile(filename)
%   channelVoltages = read_ecgfile(filename, names)
% Where:
%   channelVoltages - voltages
%   channelNames - names
%   filename - the Carto3 .txt file
%   names - is optional and is the desired channel names which helps speed

% Author: Nick Linton (2013) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------


% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

maxNVoltages = 2500;
chGain = 0.003;

filename = varargin{1};


fid = fopen(filename, 'r');
if fid == (-1)
    error(['READ_ECGFILE: Could not read the file: "' filename '"']);
end
try
    line1 = fgetl(fid);
    line2 = fgetl(fid);
    line3 = fgetl(fid);
    libChars = {'m1' 'wc'};
    if ~strcmpi(line3(1:2),libChars)
        %then we may have another version of the ecgfile
        %try the next line
        line3 = fgetl(fid);
        if ~strcmpi(line3(1:2),libChars)
            error(['READ_ECGFILE: Could not read the file: "' filename '"'])
        end
    end
    vData = fread(fid,'*char')';
    fclose(fid);
catch err
    fclose(fid);
    rethrow(err)
end
    


%line 1
spaces = isspace(line1);
line1(spaces) = [];
line1 = lower(line1);
match = strfind(line1, 'ecg_export_4.0');
if isempty(match)
    warning('READ_ECGFILE: The version number in the txt file is unexpected.') %#ok<*WNTAG>
end

%line 2
spaces = isspace(line2);
line2(spaces) = [];
line2 = lower(line2);
match = strfind(line2, 'rawecgtomv(gain)=0.003000');
if isempty(match)
    warning('READ_ECGFILE: Unexpected statement about gain.') %#ok<*WNTAG>
end

%line 3
temp = strwordpositions(line3);
nChannels = length(temp);
names = cell(1,nChannels);
temp = [temp (length(line3)+1)];
for iWord = 1:nChannels
    names{iWord} = strtrim(line3(temp(iWord):(temp(iWord+1)-1)));
end

tempChName = {};
indTempChName = 0;
i = 1;
while i<=nChannels
   
    tempStr = names{i};
    if strcmpi(tempStr(end), ')')
        tempChName{indTempChName+1} = names{i};
        indTempChName = indTempChName+1;
        i = i+1;
    else
        % look onwards from i to find a string part which ends in a ')'
        j = 0;
        nextStrPart = names{i+j};
        concatStr = names{i};
        while ~strcmpi(nextStrPart(end), ')')
            j = j + 1;
            nextStrPart = names{i+j};
            concatStr = [concatStr ' ' nextStrPart];
        end
        tempChName{indTempChName+1} = concatStr;
        %tempChName{indTempChName+1} = [names{i} ' ' names{i+1} ' ' names{i+2}];
        indTempChName = indTempChName+1;
        i = i+j+1;
    end
end
names = tempChName;
nChannels = numel(tempChName);

if nargin == 1 && nargout == 1
    varargin{1} = names; %#ok<NASGU>
    return
end

% the rest is the voltages
voltages = sscanf(vData, '%d');
if numel(voltages) ~= nChannels*2500
    error('READ_ECGFILE: Unexpected size of voltage data.')
end
voltages = chGain * reshape(voltages, nChannels, 2500)';

if nargin == 2
    desiredNames = varargin{2};
    if ischar(desiredNames)
        desiredNames = {desiredNames};
    end
    nChannels = length(desiredNames);
    k = zeros(length(desiredNames),1);
    for iCh = 1:nChannels
        k(iCh) = find(strcmpi(desiredNames{iCh}, names),1,'first');
    end
    
    % k now has the channels that we want
    voltages = voltages(:,k);
end


if nargout == 1
    varargout{1} = voltages;
elseif nargout == 2
    varargout{1} = names;
    varargout{2} = voltages;
else
    error('RED_ECGFILE: wrong number of output arguments.')
end

