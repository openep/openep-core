function success = loadensitex_prechecks(fData, data_type)
%loadensitex_prechecks Conducts checks on file headers
%   Function pulled from 'loadprecision_dxldata.m', author Steven Williams
%
% Modifications -
%   Phil Gemmell (2020): Refactored and updated
%   Steven Williams (2022): converted for EnsiteX

success = false;

% Check Export Data Element "Export Data Element : NAME"
tokens = regexp(fData, 'Export Data Element\s*:\s*(\w*)', 'once', 'tokens');
if strcmpi(data_type, 'DxL')
    goodDataElements = {'DXLData', 'DxL' };
elseif strcmpi(data_type, 'ECG')
    goodDataElements = {'ECG_FILTERED', 'ECG_RAW' };
else
    warning('No valid data type presented')
end
if ~isempty(tokens)
    dataElement = tokens{1};
    tf = strcmp(dataElement, goodDataElements);
    if all(not(tf))
        warning('LoadPrecision:InvalidFile','Invalid Export Data Element');
        return
    end
end

% Check file version
indNewLine = find(fData==newline, 1);
firstLine = fData(1:(indNewLine-1));

ind = regexp(firstLine,...
              'Export\s*File\s*Version\s*:\s*10\.0R',...
              'once');

if isempty(ind)
    warning('LoadPrecision:InvalidFile',...
            'LOADPRECISION_DXLDATA: Invalid File Revision Number');
    return
end

success = true;
end