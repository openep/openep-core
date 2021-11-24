function varargout = read_ecgfile_v4_header(varargin)
% READ_ECGFILE_V4_HEADER loads the header from an ECG file.
% Usage:
%   [electrodename_bip] = read_ecgfile_v4_header(varargin)
%   [electrodename_bip electrodename_uni] = read_ecgfile_v4_header(varargin)
%   [electrodename_bip electrodename_uni electrodename_ref] = read_ecgfile_v4_header(varargin)
% Where:
%   electrodename_bip   is the name of the electrode pair collecting the bipolar mapping point
%   electrodename_uni   is the name of the electrode collecting the unipolar mapping point
%   electrodename_ref   is the name of the reference electrode for the mapping point
%
% Author: Steven Williams (2020) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - 
%
% Info on Code Testing:
%  						 % ---------------------
%                        % test code
%                        % ---------------------
%
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

filename = varargin{1};

fid = fopen(filename, 'r');
if fid == (-1)
    error(['READ_ECGFILE_V4_HEADER: Could not read the file: "' filename '"']);
end
try
    line1 = fgetl(fid);
    line2 = fgetl(fid);
    line3 = fgetl(fid);
    if ~strcmpi(line3(1:2),'Un')

        line3 = fgetl(fid);
        if ~strcmpi(line3(1:2),'Un')
            error(['READ_ECGFILE_V4_HEADER: Could not read the file: "' filename '"'])
        end
    end
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
    warning('READ_ECGFILE_V4_HEADER: The version number in the txt file is unexpected.') %#ok<*WNTAG>
end

%line 2
spaces = isspace(line2);
line2(spaces) = [];
line2 = lower(line2);
match = strfind(line2, 'rawecgtomv(gain)=0.003000');
if isempty(match)
    warning('READ_ECGFILE_V4_HEADER: Unexpected statement about gain.') %#ok<*WNTAG>
end

%line 3 should contain the information about the mapping channels
% tokenize this line
tokens = regexp(line3, '\w*', 'match');

% % extract the unipolar; bipolar and reference electrode names
% % NOTE: THIS ONLY WORKS FOR 20A_ AND 20_B MAPPING POINTS
% uniname{1} = [tokens{4} '('];
% uniname{2} = [tokens{4}(1:4) num2str(str2double(tokens{4}(5:end))+1) '('];
% bipname = [tokens{8} '-' tokens{9}];
% refname = [tokens{12} '-' tokens{13}];

% new version
iCh = strcmpi(tokens, 'channel');
iBi = strcmpi(tokens, 'bipolar');
iRf = strcmpi(tokens, 'reference');

% Uni Name 1 is everything between 1st occurence of CHANNEL and 1st occurence of BIPOLAR
% Uni Name 2 is one greater than Uni Name 1
temp = find(iCh); iChFirst = temp(1);
temp = find(iBi); iBiFirst = temp(1);
iNeeded = iChFirst+1 : 1 : iBiFirst-1;
uniname{1} = local_concatenateName(tokens, iNeeded, ' ');
% make the second unipolar name
lengthOfFinalNumber = 2;
uni2 = '';
uniname1Cell = num2cell(uniname{1}(end-1:end));
for i = 1:numel(uniname1Cell) % check if each character is a digit
    if isnan(str2double(uniname1Cell{i}))
        lengthOfFinalNumber = lengthOfFinalNumber -1;
        continue
    else
        uni2 = [uni2 uniname1Cell{i}];
    end
end
uniname{2} = [uniname{1}(1:end-lengthOfFinalNumber) num2str(str2double(uni2)+1)];
% add on the trailing brackets
uniname{1} = [uniname{1} '('];
uniname{2} = [uniname{2} '('];

% Bipolar name is everything between 2nd occurence of CHANNEL and first occurence of REFERENCE
temp = find(iCh); iChSecond = temp(2);
temp = find(iRf); iRfFirst = temp(1);
iNeeded = iChSecond+1 : 1 : iRfFirst-1;
if strstartcmpi('MC', tokens(iNeeded(1)))
    separator = ' ';
else
    separator = '-';
end
bipname = local_concatenateName(tokens, iNeeded, separator);

% Reference name is everything after the final occurence of CHANNEL
temp = find(iCh); iChLast = temp(end);
iNeeded = iChLast+1 : 1 : numel(tokens);
refname = local_concatenateName(tokens, iNeeded, '-');

if nargout == 1
    varargout{1} = bipname;
elseif nargout == 2
    varargout{1} = bipname;
    varargout{2} = uniname;
elseif nargout == 3
    varargout{1} = bipname;
    varargout{2} = uniname;
    varargout{3} = refname;
else
    error('READ_ECGFILE_V4_HEADER: wrong number of output arguments.')
end

    function name = local_concatenateName(tokens, iNeeded, separator)
        name = '';
        for j = 1:numel(iNeeded)
            name = [name separator tokens{iNeeded(j)}];
        end
        name(1) = [];
    end

end

