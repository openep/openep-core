function S = removeEmptyFields(S)
% REMOVEEMPTYFIELDS Remove empty fields from a structure
%
% Usage:
%   S = removeEmptyFields(S)
% Where:
%   S  - the input/output structure
%
% MYFUNCTION accepts the following parameter-value pairs
%   'param1'     {value1}|vallue2
%
% MYFUNCTION Detailed description goes here
%
% Author: Steven Williams (2016)
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

% Only operate on structures
if ~isstruct(S)
    return
end

% First recurse into sub-structs
for k = 1:numel(S)
    fields = fieldnames(S(k));

    for i = 1:numel(fields)
        fname = fields{i};
        value = S(k).(fname);

        if isstruct(value)
            S(k).(fname) = removeEmptyFields(value);
        end
    end
end

% Now determine which fields are empty across ALL elements
fields = fieldnames(S);
remove = false(size(fields));

for i = 1:numel(fields)
    fname = fields{i};

    % Field is removable only if empty in every element
    remove(i) = all(arrayfun(@(x) isempty(x.(fname)), S));
end

% Remove fields in one operation (safe)
if any(remove)
    S = rmfield(S, fields(remove));
end
end