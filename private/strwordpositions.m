function varargout = strwordpositions(str)
% STRWORDSTART Finds the start of each word in a string.
% Usage:
%   iStart = strwordpositions(str)
%   [iStart, iEnd] = strwordpositions(str)
% Where:
%   iStart is a vector, or empty
%   str is a string or cell array of strings
%   iStart, iEnd are arrays or cell arrays giving the start and end of each
%                                                                     word.
%
% STRWORDSTART finds the start of each word. Where, a word is defined as
% anything that is a character or group of characters with a space before
% and after, (or the start/end of the array). The end of each word can also
% be returned.


% Author: Nick Linton (2009)
% Modifications - 

% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------

    x = isspace(str);
    if iscellstr(str)
        varargout{1} = cell(size(str));
        if nargout == 2
            varargout{2} = cell(size(str));
        end
        for i = 1:numel(str)
            if nargout == 2
                [varargout{1}{i} varargout{2}{i}] = strwordpositions(str{i});
            else
                varargout{1}{i} = strwordpositions(str{i});
            end
        end
        return
    elseif ~ischar(str);
        error('STRWORDSTART: str must be a string or cell array of strings.');
    end


    if nargout == 2
        % the end of a work is where x goes from 0 to 1
        y = ~x(1:(end-1)) & x(2:end);
        y = find(y);

        % is str(end) the end of the last word?
        if ~isspace(str(end))
            y = [y length(str)];
        end
        varargout{2} = y;
    end

    % the start of a word is where x goes from 1 to 0.
    x = x(1:(end-1)) & ~x(2:end);
    x = 1 + find(x);

    % is the str(1) the start of a word?
    if ~isspace(str(1))
        x = [1 x];
    end
    varargout{1} = x;
