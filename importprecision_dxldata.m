function data = importprecision_dxldata(varargin)
% IMPORTPRECISION_DXLDATA Loads all data from a given DxL case
%{
Imports all the ECG data from a Precision case folder, using the
loadprecision_dxldata.m file. Based in part on 'importprecision.m' file,
potentially looking to integrate later

Parameters
----------
direc : str
    Base directory in which to search for DxL data files

Flags
-----
'filematch' : cell of str
    Cell array of strings that are required to be matched in some part of
    the pathname/filename. Default={'.csv'}
'fileexclude' : cell of str
    Cell array of strings that are to be used to exclude matches, applied
    after the filematch criteria have been satisfied.
    Default={'.tif', '.fig', '.mp4', '.mpg', '.jpg', '.xml', '.png', '.surf', '.pts'}
'recursive' : bool
    true or false, whether or not the search for DxL files should be
    recursive from the base directory. Default=false
'combine' : bool
    true or false, whether to combine the individual DxL data files into
    one output for a given study. Default=true
'post_process' : bool
    true or false, whether to conduct the following post-proccessing steps:
        Zero rovLAT by min(refLAT). Default=true
'warning_files' : bool
    Whether to display warnings about being unable to load particular files
    or not. Default=false
'warning_bipole' : bool
    Whether to display warnings about determining the bipole/unipole data.
    Default=true
'warning_data' : bool
    Whether to display warnings about combining data effectively.
    Default=true

Returns
-------
data : cell of struct
    Cell array of structures, each structure containing the data for an
    individual DxL file or several combined DxL files, depending on flags.
    Contains:
        - folder
        - num_files
        - filenameList
        - study
        - sampleFreq
        - mapId
        - fileIndices
        - expStartTime
        - expEndTime
        - expStartTimeAbs
        - expEndTimeAbs
        - rovtrace
        - spare1trace
        - spare2trace
        - rovtrace_pts
        - bipole
        - ptnumber
        - rovingx
        - rovingy
        - rovingz
        - surfPtx
        - surfPty
        - surfPtz
        - rovLAT
        - peak2peak
        - peakneg
        - endtime
        - CFEmean
        - CFEstddev

Revision History
----------------
File created: Philip Gemmell (2020-02-11)
Added bipole/unipole flag: Philip Gemmell (2020-04-22)

%}

%% Check libraries have been added
search_paths = {'lib/xml_io_tools', 'lib', '../../igb_readWrite', 'panel-2.14'};
for i_path = 1:length(search_paths)
    if ~contains(path, search_paths{i_path})
        addpath(search_paths{i_path})
    end
end

%% Parse input parameters
persistent caseDirec
if isempty(caseDirec)
    caseDirec = get_homedir();
end

p = inputParser;
p.addOptional('direc',caseDirec, @(x) all(isfolder(x)));
p.addParameter('filematch',...
               {'.csv'},...
               @(x) validateattributes(x, {'string', 'cell'}, {'vector'}));
p.addParameter('fileexclude',...
               {'.tif', '.fig', '.mp4', '.mpg', '.jpg', '.xml', '.png', '.surf', '.pts'},...
               @(x) validateattributes(x, {'string', 'cell'}, {'vector'}));
p.addParameter('recursive',...
               true,...
               @islogical);
p.addParameter('combine',...
               true,...
               @islogical);
p.addParameter('post_process',...
               true,...
               @islogical);
p.addParameter('warning_files',...
               false,...
               @islogical);
p.addParameter('warning_bipole',...
               true,...
               @islogical);
p.addParameter('warning_data',...
               true,...
               @islogical);

p.addParameter('format',...
               'default',...
               @(x) validateattributes(x, {'string'} ));
p.parse(varargin{:})

if p.Results.warning_files
    warning('on', 'LoadPrecision:InvalidFile');
    warning('on', 'importprecision_dxldata:InvalidFile');
else
    warning('off', 'LoadPrecision:InvalidFile');
    warning('off', 'importprecision_dxldata:InvalidFile');
end

if p.Results.warning_bipole
    warning('on', 'importprecision_dxldata:bipolePosition');
else
    warning('off', 'importprecision_dxldata:bipolePosition');
end

if p.Results.warning_data
    warning('on', 'importprecision_dxldata:dataEntry');
else
    warning('off', 'importprecision_dxldata:dataEntry');
end

% Get working directory from inputs, if provided
if any(strcmp('direc', p.UsingDefaults))
        direc = uigetdir(p.Results.direc, 'Select the folder for the Precision Case');
    if direc == 0
        data = {};
        return
    else
        caseDirec = direc;
    end
else
    if ~iscell(p.Results.direc)
        caseDirec = {p.Results.direc};
    else
        caseDirec = p.Results.direc;
    end
end

% Adapt working directory to be recursive, if required
if p.Results.recursive
    for i_dir = 1:length(caseDirec)
        if strcmpi(caseDirec{i_dir}(end),'/')
            caseDirec{i_dir} = [caseDirec{i_dir}, '**'];
        else
            caseDirec{i_dir} = [caseDirec{i_dir}, '/**'];
        end
    end
end

[fileList, fullFileList] = get_filelist(caseDirec, p.Results.filematch,...
    p.Results.fileexclude);

if length(fileList) >= 500
    fprintf(1, "Many files found - maybe take a look?\n")
    keyboard
end

% data = cell(numel(fullFileList), 1);
data = struct;
data_bool = true(numel(fullFileList), 1);

%% Process data

% Set-up progress bar
if usejava('desktop')
    hBar = waitbar(0, 'Reading data files');
    cleanupWaitBar = onCleanup(@()close(hBar));
    set(findall(hBar, 'type', 'text'), 'Interpreter', 'none')
else
    reverseStr = '';
    fprintf('Percent done: ');
end

% Disable warnings about invalid file - we will monitor
oldWarningState = warning('query', 'loadprecision_dxldata:InvalidFile');
warning('off', 'loadprecision_dxldata:InvalidFile');
cleanupWarning = onCleanup(@()warning(oldWarningState));

for i_file = 1:length(fullFileList)
    % Extract sum total of data (while flagging for removal those entries
    % that don't provide any data)
    try
        [info, pts, egm] = loadprecision_dxldata(fullFileList{i_file});
    catch
        fprintf(1,"Can't read %s\n", fullFileList{i_file})
        continue
    end
    if isempty(info)
        warning('importprecision_dxldata:InvalidFile', ...
            ['loadprecision_dxldata: ', fileList{i_file}, ' was not loaded.'])
        data_bool(i_file) = false;
        continue
    end
    
    % Reformat data to save only ECG output
    info_fieldnames = {'study', 'sampleFreq', 'mapId', 'fileIndices', ...
        'startTime', 'endTime', 'startTimeAbs', 'endTimeAbs'};
    info_fieldnames_new = {'study', 'sampleFreq', 'mapId', 'fileIndices', ...
        'expStartTime', 'expEndTime', 'expStartTimeAbs', 'expEndTimeAbs'};
    data(i_file).filename = fullFileList{i_file};
    for iFieldname = 1:length(info_fieldnames)
        data(i_file).(info_fieldnames_new{iFieldname}) = info.(info_fieldnames{iFieldname});
    end
    
    % Save all EGM data (rovtrace, reftrace, spare1trace, spare2trace,
    % spare3trace)
    egm_elements = fieldnames(egm);
    for iEgm = 1:length(egm_elements)
        if ~strcmpi(egm_elements{iEgm}, 'reftrace')
            data(i_file).(egm_elements{iEgm}) = egm.(egm_elements{iEgm});
        end
    end
    
    % Save bipole/unipole data, with checks to make sure data all uniform.
    % Assume that all bipole data are in form "DD20   4-5", and all unipole
    % data are in form "DD20 + 12". Where such data are not recorded,
    % confirm that this is because no data recorded at all for rovtrace
    % (data may still be recorded for the ECG in spare1trace, etc., hence
    % not removing the data)
    data(i_file).rovtrace_pts = {pts(:).rovtrace};
    bipole_check = regexp(data(i_file).rovtrace_pts(:),'-');
    bipole_check = ~cellfun(@isempty, bipole_check);
    unipole_check = regexp(data(i_file).rovtrace_pts(:),'+');
    unipole_check = ~cellfun(@isempty, unipole_check);
    if all(bipole_check)
        assert(~all(unipole_check), 'Unipole signals detected in bipole');
        data(i_file).bipole = true;
    elseif all(unipole_check)
        assert(~all(bipole_check), 'Bipole signals detected in unipole');
        data(i_file).bipole = false;
    elseif any(bipole_check)
        assert(sum(unipole_check)==0, 'Both unipole and bipole signals detected')
        i_noBipoleData = find(~bipole_check);
        for i_pts = 1:length(i_noBipoleData)
            assert(all(~data(i_file).rovtrace(:,i_noBipoleData(i_pts))),...
                'Data recorded with no bipole position data.')
        end
        data(i_file).bipole = true;
        warning('importprecision_dxldata:bipolePosition',...
            ['Not all bipole positions known in ', fullFileList{i_file}])
    elseif any(unipole_check)
        assert(sum(bipole_check)==0, 'Both unipole and bipole signals detected')
        data(i_file).bipole = false;
        warning('importprecision_dxldata:bipolePosition',...
            ['Not all unipole positions known in ', fullFileList{i_file}]) 
    end
    
    % Save remaining potentially useful data
    pts_fieldnames = {'ptnumber', 'rovingx', 'rovingy', 'rovingz',...
        'surfPtx', 'surfPty', 'surfPtz',...
        'rovLAT', 'refLAT', 'peak2peak', 'peakneg', 'endtime', 'CFEmean', 'CFEstddev'};
    for iFieldname = 1:length(pts_fieldnames)
        data(i_file).(pts_fieldnames{iFieldname}) = [pts(:).(pts_fieldnames{iFieldname})];
    end
    
    if usejava('desktop')
        waitbar_msg = ['Reading data files (', num2str(i_file), '/',...
            num2str(length(fullFileList)), ')'];
        waitbar(i_file/length(fullFileList), hBar, waitbar_msg);
    else
        percentDone = 100 * i_file / length(fullFileList);
        msg = sprintf('%3.1f', percentDone);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
end
data = data(data_bool);

if isempty(data)
    fprintf(1, "No data recovered that matches requirements.\n")
    return
end
    
% Concatenate all those DxL files within a single folder
if p.Results.combine
    data = combine_data(data);
end

% Conduct requested post-processing
if p.Results.post_process
    data = post_process(data);
end

end

function data_unique = combine_data(data)

% Extract filenames and associated folder
n_data = length(data);
filelist = cell(n_data, 1);
folderlist = cell(n_data, 1);
% filelist = struct;
% folderlist = struct;
for i=1:n_data
    filelist{i} = data(i).filename;
    folder_limiter = find(filelist{i}==filesep);
    folder_limiter = folder_limiter(end);
    folderlist{i} = filelist{i}(1:folder_limiter);
end

% Confirm expected file numbers for each folder, and determine correct
% order in which to concatenate files
folderlist_unique = unique(folderlist);
n_unique = length(folderlist_unique);
data_unique = struct;
for i=1:n_unique
    data_unique(i).folder = folderlist_unique{i};
    data_unique(i).num_files = sum(contains({data.filename}, folderlist_unique{i}));
    i_match = find(contains({data.filename}, folderlist_unique{i}));
    for j=1:length(i_match)
        try
        assert(data(i_match(j)).fileIndices(2) == data_unique(i).num_files,...
            'File number mismatch!');
        catch
            keyboard
        end
    end
end
% Remove redundant section of data.fileIndices now it's been confirmed
for i=1:n_data
    data(i).fileIndices = data(i).fileIndices(1);
end
file_order = zeros(n_data,1);
i_file = 1;
for i=1:n_unique
    for j=1:data_unique(i).num_files
        file_order(i_file) = find(contains({data.filename}, folderlist_unique{i})...
            & [data.fileIndices]==j);
        i_file = i_file+1;
    end
end

% Combine data if same source folder
fieldnames_nocombine = {'filename', 'study', 'sampleFreq', 'mapId',...
    'bipole', 'fileIndices', 'expStartTime', 'expEndTime', 'expStartTimeAbs',...
    'expEndTimeAbs'};
fieldnames_nocompare = {'filename', 'fileIndices'};
folder_processed = zeros(n_unique, 1);
% dxlorder_flag = true;
% Loop over all available data
for i_order=1:n_data
    i_data = file_order(i_order);
    i_folder = find(strcmp(folderlist_unique, folderlist{i_data})==1);
    fieldname_list = fieldnames(data(i_data));
    % If folder has been processed already
    if folder_processed(i_folder)
        for i_fieldname=1:length(fieldname_list)
            
            if ~ismember(fieldname_list{i_fieldname}, fieldnames_nocombine)
                % Combine data when appropriate e.g. rovtrace
                data_unique(i_folder).(fieldname_list{i_fieldname}) = ...
                    [data_unique(i_folder).(fieldname_list{i_fieldname}), ...
                    data(i_data).(fieldname_list{i_fieldname})];
            else
                
                if ~ismember(fieldname_list{i_fieldname}, fieldnames_nocompare)
                    % Confirm that entries that are meant to be the same
                    % are actually the same e.g. sampleFreq
                    if ischar(data(i_data).(fieldname_list{i_fieldname}))
                        assert(strcmp(data(i_data).(fieldname_list{i_fieldname}),...
                            data_unique(i_folder).(fieldname_list{i_fieldname})),...
                            [fieldname_list{i_fieldname}, " doesn't match"]);
                    else
                        if ~isempty(data(i_data).(fieldname_list{i_fieldname}))
                            assert(data(i_data).(fieldname_list{i_fieldname}) == ...
                                data_unique(i_folder).(fieldname_list{i_fieldname}),...
                                [fieldname_list{i_fieldname}, " doesn't match"]);
                        else
                            warning('importprecision_dxldata:dataEntry',...
                                ['Empty entry for data_unique(',...
                                num2str(i_folder), ').', fieldname_list{i_fieldname}])
                        end
                    end
                else
                    
                    % Add filename to list of filenames for the folder
                    if strcmpi(fieldname_list{i_fieldname}, 'filename')
                        data_unique(i_folder).filenameList{end+1} = data(i_data).filename;
                        n_files = length(data_unique(i_folder).filenameList);
                        assert(data(i_data).fileIndices == n_files, 'Bugger')
                    end
                end
            end
        end
    else    % If folder hasn't been processed yet, start from scratch
        folder_processed(i_folder) = 1;
        for i_fieldname=1:length(fieldname_list)
            if ~strcmp(fieldname_list(i_fieldname), 'filename')
                data_unique(i_folder).(fieldname_list{i_fieldname}) = ...
                    data(i_data).(fieldname_list{i_fieldname});
            else
                data_unique(i_folder).filenameList = {data(i_data).filename};
            end
        end
        data_unique(i_folder).folder = folderlist_unique{i_folder};
        assert(data(i_data).fileIndices(1) == 1, 'First DxL file not recorded first!')
    end
end

assert(all(folder_processed), 'Not all unique folders processed!');

end

function data = post_process(data)

% Use points.endtime and points.rovLAT to calculate local AT
for i = 1:length(data)
    [n_data, ~] = size(data(i).rovtrace);
    assert(n_data ~= length(data(i).rovtrace_pts)); % Make sure we've got the right number!
    t_end = (n_data-1)/data(i).sampleFreq;
    data(i).rovLAT = t_end-(data(i).endtime-data(i).rovLAT);
end

end