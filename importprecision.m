function info = importprecision(varargin)
% IMPORTPRECISION loads the data from a Precision case.
% Usage:
%   info = importprecision()
%   info = importprecision(directory)
%   info = importprecision( ... , Name,Value ... )
% Output:
%   info - see 'format' below for details
% Inputs:
%   directory - an absolute folder path (if empty, use will be asked)
%   Name,Value pairs ...
%       'filematch' - a string that gives a PARTIAL match to the file(s) to
%           be loaded. e.g. 'bipol_RAW', 'Location'. Not case sensitive
%           and does NOT need to be a full match. This is useful if you do
%           not want to read in all files (save time + memory).
%       'format'
%           - 'raw' - info is a cell array with each cell
%                   corresponding to the loadfunction for each file that
%                   has been read (this is the default)
%           - 'default' - info is a structure with the following fields, which
%                      contain the info from the Data Element of similar
%                      name to the fieldname. Where there are n (n>1) files
%                      with the same Data Element, the field will have
%                      dimension n.
%                        .epcath_bip_raw
%                        .epcath_bip_filt
%                        .epcath_uni_raw
%                        .epcath_uni_filt
%                        .epcath_uni_comp
%                        .ecg_raw
%                        .ecg_filt
%                        .respiration
%                        .locations
%                        .channels
%                        .dxldata
%                        .modelgroups
%
%
% IMPORTPRECISION loads the wave data from a folder. All wavedata is
% loaded as if it is from a 'catheter' - see LOADVELOCITY_EGMDATA for more
% details.
%
% Author: Steven Williams, Nick Linton (2017)
%
% Info on Code Testing:
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

    persistent caseDirec
    if isempty(caseDirec)
        caseDirec = local_homedirec(); 
    end
    
    p = inputParser;
    p.addParameter('direc', caseDirec, @(x) isfolder(x));
    p.addParameter('filematch', {}, @(x) validateattributes(x,{'char','string','cell'}, {'vector'}));
    p.addParameter('format', 'default', @(x) validateattributes(x,{'string'} ));
    p.parse(varargin{:})
    
    if any(strcmp('direc',p.UsingDefaults))
            direc = uigetdir(p.Results.direc,'Select the folder for the Precision Case');
        if direc == 0
            info = {}; return
        else
            caseDirec = direc;
        end
    else
        caseDirec = p.Results.direc;
    end
    
    [fileList, fullFileList] = local_get_filelist(p.Results.filematch, caseDirec);

    functionList = {     @loadprecision_wavefile ...
                        ,@loadprecision_modelgroups ...
                        ,@loadprecision_dxldata ...
                        ...%,@loadprecision_shadows ...
                        ...%,@loadprecision_channels ...
                    };

    info = cell(numel(fullFileList),1);
    
    hBar = waitbar(0,'Reading data files'); 
    cleanupWaitBar = onCleanup(@()close(hBar));
    set(findall(hBar, 'type', 'text'), 'Interpreter', 'none')
    
    %disable warnings about invalid file - we will monitor
    oldWarningState = warning('query', 'LoadPrecision:InvalidFile');
    warning('off','LoadPrecision:InvalidFile');
    cleanupWarning = onCleanup(@()warning(oldWarningState));
    
    for iFile = 1:length( fullFileList )
        for iFun = 1:numel(functionList)
            temp = functionList{iFun}(fullFileList{iFile});
            if ~isempty(temp)
                info{iFile} = temp;
                break %no need to try other functions as we have succeeded
            elseif iFun == numel(functionList) %we have not loaded despite trying all functions
                warning(['IMPORTPRECISION: ' fileList{iFile} ' was not loaded.'])
            end
        end
        waitbar(iFile/length(fullFileList), hBar);
    end
    
    if isempty(info)
        return
    end
    
    info = local_reformat_precision_info(info, p.Results.format);
    info.directory = caseDirec;
    
end



function [fileList, fullFileList] = local_get_filelist(desiredFileNames, caseDirec)
    %make a list of all the '.txt. files in caseDir
    fullFileList = [];
    fileList = [];
    d = dir(caseDirec);
    if isa(desiredFileNames,'char')
        desiredFileNames = {desiredFileNames};
    end
    for iGTF = 1:length(d)
        if strcmp(d(iGTF).name,'.') || strcmp(d(iGTF).name,'..') || d(iGTF).isdir
            %do nothing
        else
            %check we have a match with listed filenames
            isToAddFile = true;
            for i = 1:numel(desiredFileNames)
                isToAddFile = false;
                startIndex = regexp(lower(d(iGTF).name), lower(desiredFileNames{i}), 'once');
                if ~isempty(startIndex)
                    isToAddFile = true;
                    break
                end
            end
            if isToAddFile
                fullFileList{end+1} = [caseDirec filesep() d(iGTF).name];        %#ok<AGROW>
                fileList{end+1} = d(iGTF).name;                                  %#ok<AGROW>
            end
        end
    end
end


function newInfo = local_reformat_precision_info(info, format)
% LOCAL REFORMAT_PRECISION_INFO reformats data from IMPORTPRECISION see below.
% Usage:
%   newInfo = reformat_precision_info(info, 'format')
% Where:
%   info is the data returned from importprecision
%   'format' is an identifier ...
%       'raw'
%       'default'
%       'SW01' - Steve Williams v01
%       'NL01' - Nick Linton v01
% Put the rest of the help above under help for IMPORTPRECISION
%
% REFORMAT_PRECISION_INFO is necessary to maintain backwards compatibility
% due to the rapidly changing nature of the SJM exports!
%
% Author: Nick Linton (2017)
%
% Info on Code Testing:
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

    %check that all info is from the same study
    
    nInfo = numel(info);
    local_checkstudy(info);
    
    switch lower(format)
        case 'raw'
            newInfo = info;
        case 'default'
            newInfo = local_default(info);
        case 'sw01'
            newInfo = local_stevewilliams_01(info);
        case 'nl01'
            newInfo = local_nicklinton_01(info);
        otherwise
            error('REFORMAT_PRECISION_INFO: format not found.')
    end

end

function newInfo = local_default(info)
% Create newInfo according to the Data Export Element. As the names of
% the Data Export Element have subtly changed some standardization is
% necessary.    
    newInfo = {};
    
    % translator =   {
    %            {'myNameA', 'precisionName1A' , 'precisionName2A' , etc} ;
    %            {'myNameB', 'precisionName1B' , 'precisionName2B' , etc} ;
    %               etc};
    
    translator = { ...
        {'epcath_bip_raw' , 'EP_Catheter_Bipolar_Raw' , 'EPcathBIObipol_RAW'};...
        {'epcath_bip_filt' , 'EP_Catheter_Bipolar_Filtered' , 'EPcathBIObipol_FILTERED'};...
        {'epcath_uni_raw' , 'EP_Catheter_Unipolar_Raw' , 'EPcathBIO_RAW'};...
        {'epcath_uni_filt' , 'EP_Catheter_Unipolar_Filtered' , 'EPcathBIO_FILTERED'};...
        {'epcath_uni_comp' , 'EPcathBIO_COMPUTED'};...
        {'ecg_raw' , 'ECG_RAW'}; ...
        {'ecg_filt' , 'ECG_FILTERED'}; ...
        {'respiration' , 'Respiration'}; ...
        {'locations' , 'Electrode_Locations' , 'Locations'}; ...
        {'channels' , 'Channels'}; ...
        {'dxldata' , '?????'}; ...
        {'modelgroups' , 'modelgroups'}; ...
                };
    
    isRead = false(numel(info),1);
    for iTr=1:numel(translator)
        isMatch = false(numel(info), 1);
        for iIn=1:numel(info)
            if isempty(info{iIn})
                isRead(iIn) = true;
            elseif isfield(info{iIn}, 'dataElement') 
                if any(contains(info{iIn}.dataElement,translator{iTr},'IgnoreCase',true))
                    isMatch(iIn) = true;
                end
            end
        end
        ind = find(isMatch);
        if numel(ind)>0
            temp = info{ind(1)};
            for i=2:numel(ind)
                temp(i)=info{ind(i)}; %concatenate the structures if more than one of same type
            end
            newInfo = setfield(newInfo,translator{iTr}{1},temp); %#ok<SFLD>
            isRead(isMatch)=true;
        end
    end
    if any(~isRead)
        beep()
        warning('REFORMAT_PRECISION_INFO: some info was not reformatted.');
    end
end

function newInfo = local_stevewilliams_01(info)
    newInfo=info;
    %     info.egmdata = egmdata; %loadprecision_egmdata.m
    %     info.modelgroups = modelgroups; %loadprecision_modelgroups.m
    %     info.shadows = shadowsdata; %loadprecision_shadows.m
    %     info.egmlocations = egmlocations; %loadprecision_electrodepositions.m
    %     info.enguidesettings = enguidesettings; %loadprecision_channels.m
    %     info.dxldata = dxldata; %loadprecision_dxldata.m        case 'nl01'
end

function newInfo = local_nicklinton_01(info)
    % Create newInfo according to the Data Export Element. As the names of
    % the Data Export Element have subtly changed some standardization is
    % necessary.
    newInfo = info;
end

function local_checkstudy(info)
    lastStudy = '';
    for i = 1:numel(info)
        if ~isempty(info{i}) && isfield(info,'study') && ~isempty(info{i}.study)
            if isempty(lastStudy)
                lastStudy = info{i}.study;
            else
                if ~strcmp(info{i}.study,lastStudy)
                    error('REFORMAT_PRECISION_INFO: files are from different studies')
                end
            end
        end
    end
end

function hd = local_homedirec()
%HOMEDIREC returns the user's home directory.

    if ispc
        hd = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    else
        hd = getenv('HOME');
    end
end