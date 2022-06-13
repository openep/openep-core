function [fileList, fullFileList] = get_filelist(caseDirec,...
    desiredFileNames, excludedFileNames)
%GET_FILELIST Make a list of all the '.txt. files in caseDir

%{
Parameters
----------
caseDirec : cell of char
    Folder in which to conduct search
desiredFileNames : cell of str
    A string that gives a PARTIAL match to the file(s) to be loaded. e.g. 
    'bipol_RAW' or 'Location'. Not case sensitive and does NOT need to be a
    full match. This is useful if you do not want to read in all files
    (save time + memory).
    NB: ** indicates doing a recursive search

Returns
-------
fileList : cell of char
    List of all matching files, without folder prefix, e.g. DxL1.csv
fullFileList : cell of char
    List of all matching files, with folder prefix, e.g. /home/dir/DxL1.csv

Revision History
----------------
2020-02-11: Philip Gemmell
    Extracted from importprecision.m. Restructured to include preallocation

%}

% Create initial list of files, which will be removed from subject to
% requirements
if ~iscell(caseDirec)
    caseDirec = {caseDirec};
end
n_dir = length(caseDirec);
% d = cell(n_dir, 1);
fullFileList = cell(n_dir, 1);
fileList = cell(n_dir, 1);
for i_dir = 1:n_dir
    d = dir(caseDirec{i_dir});
    fullFileList{i_dir} = cell(length(d), 1);
    fileList{i_dir} = cell(length(d), 1);
    match_bool = false(length(d), 1);

    for i_file = 1:length(d)
        fullFileList{i_dir}{i_file} = [d(i_file).folder, filesep(), d(i_file).name];
        fileList{i_dir}{i_file} = d(i_file).name;

        if strcmp(d(i_file).name, '.') || strcmp(d(i_file).name, '..') || d(i_file).isdir
            % Do nothing
        else
            % Check we have a match with listed filenames (if any provided)
            match_bool(i_file) = true;
            for i_filematch = 1:numel(desiredFileNames)
    %             match_bool(i_file) = false;
                startIndex = regexpi(fullFileList{i_dir}{i_file}, desiredFileNames{i_filematch}, 'once');
                if isempty(startIndex)
                    match_bool(i_file) = false;
                    break
                end
    %             if ~isempty(startIndex)
    %                 match_bool(i_file) = true;
    %                 break
    %             end
            end

            % Check to exclude any filenames (if required)
            for i_fileexclude = 1:numel(excludedFileNames)
                startIndex = regexpi(fullFileList{i_dir}{i_file}, excludedFileNames{i_fileexclude}, 'once');
                if ~isempty(startIndex)
                    match_bool(i_file) = false;
                    break
                end
            end

        end
    end
    fullFileList{i_dir} = fullFileList{i_dir}(match_bool);
    fileList{i_dir} = fileList{i_dir}(match_bool);
end

fullFileList = vertcat(fullFileList{:});
fileList = vertcat(fileList{:});
end