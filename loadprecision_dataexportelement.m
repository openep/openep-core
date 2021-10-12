function dee = loadprecision_dataexportelement(filename)
% LOADPRECISION_DATAEXPORTELEMENT reads the data export element from a file
% Usage:
%   dee = loadprecision_dataexportelement(filename)
% Where:
%   filename is the filename
%   dee is the output

% Author: Nick Linton (2017)
% Modifications
%
% Info on Code Testing:
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------


% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

    fileID = fopen(filename, 'r');
    if fileID == (-1); error('LOADPRECISION_DATAEXPORTELEMENT: Could not open file.'); end
    cleanupFile = onCleanup(@()fclose(fileID));

    % read the start data into the memory 
    maxBytes = 1000; % enough data to cover the start of the header   
    fseek(fileID, 0, 'bof');
    if maxBytes > filebytes2end(fileID); maxBytes = filebytes2end(fileID); end
    [fData, fDataSize] = fread(fileID, maxBytes, '*char');
    fData = fData(1:fDataSize)';
    
    tokens = regexp(fData, 'Export Data Element\s*:\s*(\w*)', 'once', 'tokens');
    if ~isempty(tokens)
        dee = tokens{1};
    else
        dee = '';
    end
end
