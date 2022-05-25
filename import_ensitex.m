function [userdata, matFileFullPath] = import_ensitex(varargin)
% IMPORT_ENSITEX is used to import an EnsiteX case.
%
% Usage:
%   userdata = import_ensitex()
%   userdata = import_ensitex(directory)
%   userdata = import_ensitex( ... , Name,Value ... )
%
% Where:
%   directory - an absolute folder path (if empty, user will be asked)
%   userdata - an OpenEP data structure
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

% Parse input data
persistent saveDir homeDir
if isempty(saveDir) || ~ischar(saveDir)
    saveDir = userpath();
end
if ~isfolder(saveDir)
    saveDir = userpath();
end

if isempty(homeDir)
    homeDir = userpath()';
end
if ~isfolder(homeDir)
    homeDir = saveDir;
end

userdata = [];
hWait = [];

if nargin >= 1
    userinput = varargin{1};
else
    dialog_title = 'Select the EnsiteX Contact_Mapping_Model file, or a .mat file.';
    filterSpec = {'*.xml;*.mat', 'Appropriate files (*.xml;*.mat)' ; '*.xml','XML (*.xml)' ; '*.mat','Matlab (*.mat)' ; '*.*','All files (*.*)'};
    if ~ispc()
        uiwait(msgbox(dialog_title,'modal'))
    end
    [filename,pathname] = uigetfile(filterSpec, dialog_title, homeDir);
    userinput = fullfile(pathname, filename);
    if filename == 0
        return
    else
        homeDir = pathname;
    end
end

if strcmpi(userinput((end-3):end),'.zip')
    %then we need to unzip the folder
    error('IMPORT_ENSITEX: ZIP files are no longer supported - unzip externally to Matlab.');
    return %#ok<UNRCH>
elseif strcmpi(userinput((end-3):end),'.mat')
    s = load(userinput);
    userdata = s.userdata;
    return
else
    studyDir = fileparts(userinput);
    homeDir = studyDir;
end

% parse command line input
nStandardArgs = 1;
bipoleType = 'along';
channelRef_cli = '';
channelECG_cli = '';
saveFileName_cli = '';
verbose = true;
if nargin > nStandardArgs
    for i = nStandardArgs+1:2:nargin
        switch lower(varargin{i})
            case 'bipoletype'
                bipoleType_cli = varargin{i+1};
            case 'refchannel'
                channelRef_cli = varargin{i+1};
                if ischar(channelRef_cli)
                    channelRef_cli = {channelRef_cli};
                end
            case 'ecgchannel'
                channelECG_cli = varargin{i+1};
                if ischar(channelECG_cli)
                    channelECG_cli = {channelECG_cli};
                end
            case 'savefilename'
                saveFileName_cli = varargin{i+1};
            case 'verbose'
                verbose = varargin{i+1};
            otherwise
                error('IMPORTCARTO_MEM: Unrecognized input.')
        end
    end
end


% Create an empty OpenEP data structure
userdata = openep_createuserdata;

% General data
userdata.systemName = 'ensitex';
userdata.notes{1} = [date() ': Created'];
userdata.ensiteXFolder = studyDir;

% Load the model groups
info = loadprecision_modelgroups([studyDir filesep() 'Contact_Mapping_Model.xml']);

% assign the geometry data
if isfield(info.dxgeo, 'triangles') && isfield(info.dxgeo, 'vertices')
    TRI = info.dxgeo.triangles;
    X = info.dxgeo.vertices(:,1);
    Y = info.dxgeo.vertices(:,2);
    Z = info.dxgeo.vertices(:,3);
    userdata.surface.triRep = TriRep(TRI, X, Y, Z);
    FF = freeBoundary(userdata.surface.triRep);
    isVertexAtRim = false(size(userdata.surface.triRep.X,1),1);
    if ~isempty(FF)
        isVertexAtRim(FF(:,1)) = true;
    end
    userdata.surface.isVertexAtRim = isVertexAtRim;
end

% deal with labels, for the test case this is 0, 1 or 2
allLabels = info.dxgeo.surface_of_origin;
labels = unique(allLabels);
cMap = colormap(parula(numel(labels)));
faceColors = trFaceToVertData(userdata.surface.triRep, allLabels);
hSurf = drawMap(userdata, 'type', 'none');
colorShell(hSurf, [X Y Z], faceColors, Inf ...
    , 'showcolorbar', 'show' ...
    , 'coloraxis', [0 numel(labels)] ...
    , 'interpolation', 'off' ...
    , 'usrColorMap', cMap ...
    , 'datatype', 'labels' ...
    );

% next load the map file and the wave files into memory using loadensitex_dxldata.m
cMapDir = [studyDir filesep() 'Contact_Mapping'];
mappingFiles = nameFiles(cMapDir, 'showhiddenfiles', false, 'extension', '.csv');
for i = 1:numel(mappingFiles)
    [info, varnames, data] = loadensitex_dxldata([varargin{1} filesep() 'Contact_Mapping' filesep() mappingFiles{i}]);
    dataFile{i}.info = info; %#ok<*AGROW> 
    dataFile{i}.varnames = varnames;
    dataFile{i}.data = data;
end

% s = load('/Users/steven/Desktop/dataFiles.mat');
% dataFile = s.dataFile;

% find the map name
iContactMapFile = local_findFile('Map_CV_omni.csv');
map.name = dataFile{iContactMapFile}.info.mapName;
map.type = dataFile{iContactMapFile}.info.mapType;
map.study = dataFile{iContactMapFile}.info.study;

% then map the data into the OpenEP data format

% Deal wtih the infoMapping dictionary
infoMapping = { ...
    'electric.sampleFrequency'                'Wave_rov.csv'                              'sampleFreq' ...
    ; ...
    };
for i = 1:size(infoMapping)
    fieldNames = strsplit(infoMapping{i,1}, '.');
    thisFileName = infoMapping{i,2};
    fileInd = local_findFile(thisFileName, dataFile);
    thisFieldName = infoMapping{i,3};
    switch numel(fieldNames)
        case 1
            userdata.(fieldNames{1}) = dataFile{fileInd}.info.(thisFieldName);
        case 2
            userdata.(fieldNames{1}).(fieldNames{2}) = dataFile{fileInd}.info.(thisFieldName);
        case 3
            userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}) = dataFile{fileInd}.info.(thisFieldName);
        case 4
            userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}).(fieldNames{4}) = dataFile{fileInd}.info.(thisFieldName);
        otherwise
            error('OPENEP/IMPORT_ENSITEX: Code not yet implemented for more than 4 sub fields')
    end
end

% I think we have to decide whether to import an 'along' map or an 'across'
% map and these might need combining retrospectively for particular use
% cases, but could be considered as separate OpenEP data structures at
% import time. The default should possibly be to import to the 'along' as
% this is 'along the spline' and consistent with usual practice in EP,
% whereas 'across' may be influenced by HD grid geometry.

% Deal wtih the dataMapping dictionary
dataMapping = { ...
       'electric.tags'                        'Map_CV_omni.csv'                           'annot' ...
    ;  'electric.names'                       'Map_CV_omni.csv'                           '(Point #)' ...
    ;  'electric.electrodeNames_bip'          'Wave_bi_along.csv'                         'Trace' ...
    ;  'electric.egmX'                        'Map_CV_omni.csv'                           'roving x,roving y,roving z' ...
    ;  'electric.egm'                         'Wave_bi_along'                             'signals' ...
    ;  'electric.electrodeNames_uni'          'Wave_uni_along'                            'Trace' ...
    ;  'electric.egmUniX'                     'Map_CV_omni.csv'                           'Uni_CornerX,Uni_CornerY,Uni_CornerZ,Uni_AlongX,Uni_AlongY,Uni_AlongZ' ...
    ;  'electric.egmUni'                      'Wave_uni_corner.csv,Wave_uni_along.csv'    'signals' ...
    ;  'electric.egmRef'                      'Wave_refs.csv'                             'signals' ... % Ask, 'which signal is a good reference'
    ;  'electric.ecg'                         'Wave_refs.csv'                             'signals' ... % Ask, 'which signal is a good ECG'
    ;  'electric.annotations.woi'             'Map_CV_omni.csv'                           'left curtain (ms),right curtain (ms)' ... % this is in ms relative to samples!
    ;  'electric.annotations.referenceAnnot'  'Map_CV_omni.csv'                           'Ref Tick' ... % this is in samples!
    ;  'electric.annotations.mapAnnot'        'Wave_bi_along.csv'                         'rovTime (wave samples)' ... % this is _presumably_ in samples!
    ;  'electric.voltages.bipolar'            'Map_CV_omni.csv'                           'pp_Valong' ...
    ;  'electric.voltages.unipolar'           'Map_CV_omni.csv'                           'unipoleMaxPP' ...
    ;  'electric.egmSurfX'                    'Map_CV_omni.csv'                           'surface x,surface y,surface z' ...
    ;  'electric.barDirection'                'Map_CV_omni.csv'                           'normal x,normal y,normal z' ...
    ;  'electric.include'                     'Map_CV_omni.csv'                           'utilized' ...
    };

   % ;  'electric.impedances.time'             ''                                          '' ... % We do not seem to have impedance data
   % ;  'electric.impedances.value'            ''                                          '' ... % We do not seem to have impedance data

for i = 1:size(dataMapping)
    
    fieldNames = strsplit(dataMapping{i,1}, '.');
    thisFileName = dataMapping{i,2};
    fileInd = local_findFile(thisFileName, dataFile);

    thisFieldName = dataMapping{i,3};
    
    fieldInd = local_findField(thisFieldName, dataFile(fileInd));
    % parse the data

    switch dataMapping{i,1}
        case 'electric.tags' 
            userdata.electric.tags                          = dataFile{fileInd}.data(:,fieldInd);

        case 'electric.names'  
            userdata.electric.names                         = dataFile{fileInd}.data(:,fieldInd);

        case 'electric.electrodeNames_bip' 
            userdata.electric.electrodeNames_bip            = dataFile{fileInd}.data(:,fieldInd);

        case 'electric.egmX' 
            X = str2double(dataFile{fileInd}.data(:,fieldInd(1)));
            Y = str2double(dataFile{fileInd}.data(:,fieldInd(2)));
            Z = str2double(dataFile{fileInd}.data(:,fieldInd(3)));
            userdata.electric.egmX                          = [X Y Z];

        case 'userdata.electric.egm'
            userdata.electric.egm                           = cell2mat(dataFile{fileInd}.data(:,fieldInd));

        case 'electric.electrodeNames_uni'
            userdata.electric.electrodeNames_uni            = dataFile{fileInd}.data(:,fieldInd);

        case 'electric.egmUniX'
            X = str2double(dataFile{fileInd}.data(:,fieldInd(1)));
            Y = str2double(dataFile{fileInd}.data(:,fieldInd(2)));
            Z = str2double(dataFile{fileInd}.data(:,fieldInd(3)));
            userdata.electric.egmUniX(:,:,1)                = [X Y Z];

        case 'electric.egmUni'
            userdata.electric.egmUni(:,:,1)                 = cell2mat(dataFile{fileInd(1)}.data(:,fieldInd));
            userdata.electric.egmUni(:,:,2)                 = cell2mat(dataFile{fileInd(2)}.data(:,fieldInd));

        case 'electric.egmRef'
            % Hard coding a channel just now; but ask 'which signal is a good reference?' of the user
            userdata.electric.egmRef                        = userdata.electric.egm; % TODO only for now

        case 'electric.ecg'
            % Hard coding a channel just now; but ask 'which signal is a good reference?' of the user
            userdata.electric.ecg                           = userdata.electric.egm; % TODO only for now

        case 'electric.annotations.woi'
            startWindow = str2double(dataFile{fileInd}.data(:,fieldInd(1)));
            endWindow = str2double(dataFile{fileInd}.data(:,fieldInd(2)));
            userdata.electric.annotations.woi               = round([startWindow endWindow] ./ 1000 .* userdata.electric.sampleFrequency);

        case 'electric.annotations.referenceAnnot'
            userdata.electric.annotations.referenceAnnot    = str2double(dataFile{fileInd}.data(:,fieldInd));

        case 'electric.annotations.mapAnnot'  
            userdata.electric.annotations.mapAnnot          = round(str2double(dataFile{fileInd}.data(:,fieldInd)) ./ userdata.electric.sampleFrequency .* userdata.electric.sampleFrequency);

        case 'electric.voltages.bipolar' 
            userdata.electric.voltages.bipolar              = str2double(dataFile{fileInd}.data(:,fieldInd));

        case 'electric.voltages.unipolar' 
            userdata.electric.votlages.unipolar             = str2double(dataFile{fileInd}.data(:,fieldInd));

        case 'electric.impedances.time'
            % TODO

        case 'electric.impedances.value' 
            % TODO

        case 'electric.egmSurfX'  
            X = str2double(dataFile{fileInd}.data(:,fieldInd(1)));
            Y = str2double(dataFile{fileInd}.data(:,fieldInd(2)));
            Z = str2double(dataFile{fileInd}.data(:,fieldInd(3)));
            userdata.electric.egmSurfX                      = [X Y Z];

        case 'electric.barDirection' 
            X = str2double(dataFile{fileInd}.data(:,fieldInd(1)));
            Y = str2double(dataFile{fileInd}.data(:,fieldInd(2)));
            Z = str2double(dataFile{fileInd}.data(:,fieldInd(3)));
            userdata.electric.barDirection                  = [X Y Z];

        case 'electric.include'
            userdata.electric.include                       = str2double(dataFile{fileInd}.data(:,fieldInd));
    end

    % save the fieldnames
%     switch numel(fieldNames)
%         case 1
%             userdata.(fieldNames{1}) = dataToSave;
%         case 2
%             userdata.(fieldNames{1}).(fieldNames{2}) = dataToSave;
%         case 3
%             userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}) = dataToSave;
%         case 4
%             userdata.(fieldNames{1}).(fieldNames{2}).(fieldNames{3}).(fieldNames{4}) = dataToSave;
%         otherwise
%             error('OPENEP/IMPORT_ENSITEX: Code not yet implemented for more than 4 sub fields')
%     end
end

    function iF = local_findFile(f, d)
        % find the index into the cell array, d, of the filename f
        for iD = 1:numel(d)
            [~,n,e] = fileparts(d{iD}.info.filename);
            allFileNames{iD} = [n e]; %#ok<AGROW>
        end
        requiredFiles = strsplit(f,',');
        for iFile = 1:numel(requiredFiles)
            iF(iFile) = find(strstartcmpi(requiredFiles{iFile},allFileNames)); %#ok<AGROW> 
        end
    end

    function iC = local_findField(f, d)
        % find the column iC in d.data such that d.varname{iC} == f
        requiredFields = strsplit(f,',');
        for iDataFile = 1:numel(d)
            for iField = 1:numel(requiredFields)
                iC(iDataFile,iField) = find(strcmpi(d{iDataFile}.varnames,requiredFields{iField})); %#ok<AGROW>
            end
        end
    end

% userdata.electric.sampleFrequency = dataFile{5}.info.sampleFreq;

% userdata.electric.tags                          = dataFile{1}.data(:,26);
% userdata.electric.names                         = dataFile{1}.data(:,5);
% userdata.electric.electrodeNames_bip            = dataFile{3}.data(:,1);
% userdata.electric.egmX                          = [str2double(dataFile{1}.data(:,7)) str2double(dataFile{1}.data(:,8)) str2double(dataFile{1}.data(:,9))];
% userdata.electric.egm                           = cell2mat(dataFile{3}.data(:,6));
% userdata.electric.electrodeNames_uni            = dataFile{7}.data(:,1);
% userdata.electric.egmUniX(:,:,1)                = [str2double(dataFile{1}.data(:,60)) str2double(dataFile{1}.data(:,61)) str2double(dataFile{1}.data(:,62))];
% userdata.electric.egmUniX(:,:,2)                = [str2double(dataFile{1}.data(:,64)) str2double(dataFile{1}.data(:,65)) str2double(dataFile{1}.data(:,66))];
% userdata.electric.egmUni(:,:,1)                 = cell2mat(dataFile{8}.data(:,6));
% userdata.electric.egmUni(:,:,2)                 = cell2mat(dataFile{7}.data(:,6));
% % Hard coding a channel just now; but ask 'which signal is a good reference?' of the user
% userdata.electric.egmRef                        = userdata.electric.egm; % only for now
% % Hard coding a channel just now; but ask 'which signal is a good reference?' of the user
% userdata.electric.ecg                           = userdata.electric.egm; % only for now
% userdata.electric.annotations.woi               = round([str2double(dataFile{1}.data(:,23)) str2double(dataFile{1}.data(:,24))] ./ 1000 .* userdata.electric.sampleFrequency);
% userdata.electric.annotations.referenceAnnot    = str2double(dataFile{1}.data(:,74));
% userdata.electric.annotations.mapAnnot          = round(str2double(dataFile{3}.data(:,5)) ./ userdata.electric.sampleFrequency .* userdata.electric.sampleFrequency);
% userdata.electric.voltages.bipolar              = str2double(dataFile{1}.data(:,32));
% userdata.electric.votlages.unipolar             = str2double(dataFile{1}.data(:,58));
% userdata.electric.egmSurfX                      = [str2double(dataFile{1}.data(:,10)) str2double(dataFile{1}.data(:,11)) str2double(dataFile{1}.data(:,12))];
% userdata.electric.barDirection                  = [str2double(dataFile{1}.data(:,13)) str2double(dataFile{1}.data(:,14)) str2double(dataFile{1}.data(:,15))];
% userdata.electric.include                       = str2double(dataFile{1}.data(:,17));



% % Ablation data - TODO
% % userdata.rf.originaldata.force.time =
% % userdata.rf.originaldata.force.force =
% % userdata.rf.originaldata.force.axialangle =
% % userdata.rf.originaldata.force.lateralangle =
% % userdata.rf.originaldata.force.position =
% % userdata.rf.originaldata.ablparams.time =
% % userdata.rf.originaldata.ablparams.power =
% % userdata.rf.originaldata.ablparams.impedance =
% % userdata.rf.originaldata.ablparams.distaltemp =


    % Encourage user to save the data
    if ~isempty(saveFileName_cli)
        save(saveFileName_cli, 'userdata');
        matFileFullPath = saveFileName_cli;
    else
        defaultName = [map.study '_' map.name];
        defaultName(isspace(defaultName)) = '_';
        originalDir = cd();
        matFileFullPath = fullfile(saveDir, defaultName); %default
        cd(saveDir);
        [filename,saveDir] = uiputfile('*.mat', 'Save the userdata to disc for future rapid access?',defaultName);
        cd(originalDir);
        if filename ~= 0
            save([saveDir filename], 'userdata','-v7.3'); %needed as sometimes >2GB
            matFileFullPath = fullfile(saveDir, filename);
        end
    end


end