function [x1, f_x1] = removeDuplicateDataPoints(x, f_x)
% REMOVEDUPLICATEDATAPOINTS removes duplicate datapoints in x and f_x
%
% Usage:
%   [x1, f_x1] = removeDuplicateDataPoints(x, f_x)
% Where:
%   x       - locations
%   f_x     - data values
%   x1      - locations with duplicates removed
%   f_x1     - new data values
%
% REMOVEDUPLICATEDATAPOINTS Detects duplicate data points in x and removes
% these. The algorithm averages the corresponding values in f_x; unless
% f_x(n) == min(f_x) || f_x(n) = max(f_x) in which case the minimum or
% maximum values are respected.
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

[u,I,J] = unique(x, 'rows', 'first');
hasDuplicates = size(u,1) < size(x,1);

x1 = x;
f_x1 = f_x;

if hasDuplicates
    ixDupRows = setdiff(1:size(x,1), I);
    rowsThatHaveDuplicates = unique(J(ixDupRows));

    min_f_x = min(f_x);
    max_f_x = max(f_x);

    % we want to modify the value in f_x of the first instance of a
    % duplicated row in x; we will remove the duplicates later
    iRemove = [];
    for i = 1:numel(rowsThatHaveDuplicates)

        % find the set of values in f_X
        iTheseDuplicates = J == rowsThatHaveDuplicates(i);
        theseValues = f_x(iTheseDuplicates);

        % average them (or keep the min or max value)
        if min(theseValues) == min_f_x
            newValue = min(theseValues);
        elseif max(theseValues) == max_f_x
            newValue = max(theseValues);
        else
            newValue = mean(theseValues);
        end

        if min(theseValues) == min_f_x && max(theseValues) == max_f_x
            newValue = mean(theseValues);
            warning('OPENEP:information','OPENEP/REMOVEDUPLICATEDATAPOINTS: A duplicate is detected at both the global minimum and global maximum of f_x');
        end

        % assign this value to the first occurence
        f_x1(rowsThatHaveDuplicates(i)) = newValue;

        % save the set of numbers to remove
        duplicateLocations = find(iTheseDuplicates);
        iRemove = [iRemove; duplicateLocations(2:end)];
    end

    x1(iRemove,:) = [];
    f_x1(iRemove) = [];

end

end