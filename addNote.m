function userdata = addNote(userdata, note, varargin )
% ADDNOTE Adds notes to an OpenEP dataset
%
% Usage:
%   userdata = addNote(userdata, note, varargin)
% Where:
%   userdata  - the output
%   userdata  - the input
%   note  - note string
%
% ADDNOTE accepts the following parameter-value pairs
%   'type'     {'date'}|'label'
%
% By default, ADDNOTE adds notes with dates. Alternatively, the parameter
% 'type' can be specified to be a label, in which case a note with this
% label will be added
%
% Author: Steven Williams (2022)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 2; % UPDATE VALUE
type = 'date';
label = [];
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'type'
                type = varargin{i+1};
            case 'label'
                label = varargin{i+1};
        end
    end
end

switch type
    case 'date'
        userdata.notes{end+1} = [date ': ' note];

    case 'label'
        if isempty(label)
            error('OPENEP/ADDNOTE: Note with type LABEL specified but no LABEL provided');
        else
            userdata.notes{end+1} = [label ': ' note];
        end
end

end