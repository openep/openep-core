function tf = inserttextintotextfile(fName, line2insert, extraText)
% INSERTTEXTINTOTEXTFILE Inserts text into a .txt file.
%
% Usage:
%   tf = inserttextintotextfile(fid, line2insert, extraText)
%
% Where:
%   tf true if successful
%   fName is the file name
%   line2insert is the line at which text is inserted
%   extraText is the text to insert (a newline is added if needed)


% Author: Nick Linton (2013)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

% check that the extraText starts with a newline
if ~ischar(extraText)
    error('ISERTTEXTINTOTEXTFILE: extraText should be char')
end
if extraText(1)~= char(10)
    extraText = [char(10) extraText];
end
        
tf = false; 

fid = fopen(fName, 'rt+');
if fid<3
    error('ISERTTEXTINTOTEXTFILE: unable to open file')
end

%read the whole file into memory
frewind(fid)
data = fread(fid, inf, '*char');

%find the position where we are going to break the text.
if isnumeric(line2insert)
    if line2insert <= 1;
        line2insert = 'first';
    else
        nLine = 1;
        insertPosition = 0;
        while (nLine < line2insert)
            insertPosition = insertPosition+1;
            if insertPosition >= numel(data)
                break
            elseif data(insertPosition+1) == char(10)
                nLine = nLine+1;
            end
        end
    end
if ischar(line2insert)
    if strcmpi(line2insert, 'last')
        insertPosition = numel(data);
    elseif strcmpi(line2insert, 'first')
        insertPosition = 0;
        extraText(1) = [];
        if extraText(end) ~= char(10)
            extraText = [extraText char(10)];
        end
    else
        error('unknown command')
    end
end
        
%now write the file
fseek(fid, 0, 'bof');
fwrite(fid, data(1:(insertPosition)));
fwrite(fid, extraText);
fwrite(fid, data((insertPosition+1):end));

fclose(fid);

tf = true;
end