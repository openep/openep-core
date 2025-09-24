function [userdata, matFileFullPath] = importensitex_openep(varargin)
% IMPORTPRECISION_OPENEP provides a data structure from multiple Precision files.
% Usage:
%   userdata = importprecision_openep(userinput)
%   userdata = importprecision_openep()
%   [userdata, matFileFullPath] = ...
% Where:
%   dirName is the directory with all of the files corresponding to a map
%   userdata is a single data structure
%   matFileFullPath is the path to the .mat file, if opened or saved
%
% IMPORTPRECISION_OPENEP accepts the following parameter-value pairs:
%   'savefilename'      {''}|string
%       The full path to the location in which to save the output.
%   'type'              {'standard'}|'omnipolar'
%       Specifies whether to import a standard map, based on local
%       activation time (LAT), bipolar (BIP) and unipolar (UNI) folders, or
%       whether to important an omnipolar dataset
%
% Example of command line entry ...
%       userdata = importensitex_openep(<path to folder>, ...
%                                        'savefilename',
%                                        '~/Desktop/test.mat');
%
% userdata structure ...
%   .surface
%       .triRep         - TriRep object for the surface
%       .isVertexAtRim  - logical array indicating vertices at a 'rim'
%       .act_bip        - nVertices*2 array of activation and voltage data
%       .uni_imp_frc    - nVertices*3 array of uni voltage, impedance and contact force
%   .electric
%       .isPointLocationOnly    - logical array
%       .tags
%       .names
%       .egmX           - location of point
%       .egmSurfX       - location of surface nearest point
%       .barDirection   - normal to surface at egmSurfX
%       .egm            - bipolar electrogram
%       .egmUni         - matrix of unipolar electrograms
%       .egmUniX        - location of unipolar points
%       .egmRefNames    - names of egmRef
%       .egmRef         - electrogram of reference
%       .ecgNames       - ecg names (or other channel names)
%       .ecg            - ecg
%       .force
%           .force    - instantaneous force recording
%           .axialAngle    - axial angle
%           .lateralAngle  - lateral angle
%           .time_force - time course of force [(:,:,1)=time, (:,:,2)=force]
%           .time_axial - time course of axial angle [(:,:,1)=time, (:,:,2)=axial angle]
%           .time_lateral - time course of lateral angle [(:,:,1)=time, (:,:,2)=lateral angle]

% Author: Steven Williams (2025)
% SPDX-License-Identifier: Apache-2.0
%
% ---------------------------------------------------------------
% testing
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

%% Identify the folder containing the exported files
persistent saveDir homeDir
if isempty(saveDir) || ~ischar(saveDir) || ~isfolder(saveDir)
    saveDir = local_homedirec();
end
if isempty(homeDir) || ~isfolder(homeDir); homeDir = saveDir; end
if nargin >= 1
    userinput = varargin{1};
else
    dialog_title = 'Select the folder containing the exported EnsiteX files';
    if ~ispc()
        uiwait(msgbox(dialog_title,'modal'))
    end
    userinput = uigetdir(homeDir, dialog_title);
    if userinput == 0
        return
    else
        homeDir = userinput;
    end
end
if ~isfolder(userinput)
    error('IMPORTENSITEX_OPENEP: Only specify a folder as user input')
else
    studyDir = userinput;
    homeDir = studyDir;
end

%% Parse command line input
nStandardArgs = 1; % UPDATE VALUE
saveFileName_cli = '';
if nargin > nStandardArgs
    for i = nStandardArgs+1:2:nargin
        switch varargin{i}
            case 'savefilename'
                saveFileName_cli = varargin{i+1};
            otherwise
                error('IMPORTCARTO_MEM: Unrecognized input.')
        end
    end
end

%% Identify the relevant subfolders

% These should be the BIP folder, the LAT folder and the UNI folder
latMapDir = local_findDirectory('OpenEP_LAT', studyDir);
bipMapDir = local_findDirectory('OpenEP_BIP', studyDir);
uniMapDir = local_findDirectory('OpenEP_UNI', studyDir);
omniPPMapDir = local_findDirectory('OpenEP_OMNI_PP', studyDir);
omniLATMapDir = local_findDirectory('OpenEP_OMNI_LAT', studyDir);

% And including the set of ALL the automap folders
allSubFolds = nameFolds(studyDir);
automapDirs = allSubFolds(strstartcmpi('OpenEP_AutoMap', allSubFolds));

%% Parse the geometry and surface mapping data
% By loading the LAT XML file to get the geometry
% By loading the LAT XML file to get the LAT map (I know, unnecessary
% duplication but this gives us future flexibility)
% By loading the bipolar XML file to get the bipolar map
% By loading the unipolar XML file to get the unipolar map
data_geometry = loadprecision_modelgroups(fullfile(latMapDir{:}, 'Contact_Mapping_Model.xml'));
TRI = data_geometry.dxgeo.triangles;
X = data_geometry.dxgeo.vertices(:,1);
Y = data_geometry.dxgeo.vertices(:,2);
Z = data_geometry.dxgeo.vertices(:,3);
tr = TriRep(TRI, X, Y, Z);
t.X = tr.X;
t.Triangulation = tr.Triangulation;
normals = data_geometry.dxgeo.normals;

data_latMap = loadprecision_modelgroups(fullfile(latMapDir{:}, 'Contact_Mapping_Model.xml'));
data_bipolarVoltageMap = loadprecision_modelgroups(fullfile(bipMapDir{:}, 'Contact_Mapping_Model.xml'));
data_unipolarVoltageMap = loadprecision_modelgroups(fullfile(uniMapDir{:}, 'Contact_Mapping_Model.xml'));
act = data_latMap.dxgeo.act;
%act(data_latMap.dxgeo.map_status==2) = NaN;
bip = data_bipolarVoltageMap.dxgeo.bip;
%bip(data_bipolarVoltageMap.dxgeo.map_status==2) = NaN;
uni = data_unipolarVoltageMap.dxgeo.bip;
%uni(data_unipolarVoltageMap.dxgeo.map_status==2) = NaN;
imp = NaN(size(uni));
frc = NaN(size(uni));

%% (Optional) Parse additional voltage maps - e.g. omnipolar
% By loading the XML file

%% Parse electrogram data
% By loading the DXL files (note that these differ for unipolar and bipolar)

%load the mapping points data
[info, varnames, data] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() 'Map_LAT_bi.csv']);
mappingPoints.info = info; 
mappingPoints.varnames = varnames; 
mappingPoints.data = data;

%additionally get the measured voltages for each point
[info, varnames, data] = loadensitex_dxldata([bipMapDir{:} filesep() 'Contact_Mapping' filesep() 'Map_PP_bi.csv']);
info = info; 
varnames = varnames;
data = data;

% do some simple checks for compatibility between the PP and LAT datasets
isError = false;
if mappingPoints.info.numPoints ~= info.numPoints
    warning('IMPORTENSITEX_OPENEP: Mismatch between number of points in the voltage and activation time datasets');
    isError = true;
end
if ~strcmpi(mappingPoints.info.mapName, info.mapName)
    warning('IMPORTENSITEX_OPENEP: Mismatch between map names in the voltage and activation time datasets');
    isError = true;
end
if ~strcmpi(mappingPoints.info.study, info.study)
    warning('IMPORTENSITEX_OPENEP: Mismatch between study names in the voltage and activation time datasets');
    isError = true;
end
if isError
    error('IMPORTENSITEX_OPENEP: Error parsing data. See warnings above for hints');
end
% TODO: there are likely to be other checks we could add in here

% access the voltage data from the PP data and save along with the LAT data
ppStr = 'P-P';
ppValidStr = 'P-P valid';
mappingPoints.varnames{end+1} = ppStr;
mappingPoints.varnames{end+1} = ppValidStr;
ppData = data(:,strcmpi(varnames,ppStr));
ppValidData = data(:,strcmpi(varnames,ppValidStr));

% concatenate
mappingPoints.data = [mappingPoints.data ppData ppValidData];

%get the reference electrograms
[refInfo, refVarnames, refData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() 'Wave_refs.csv']);

%get the roving bipolar electrograms
[rovInfo, rovVarnames, rovData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() 'Wave_rov.csv']);

%get the roving unipolar distal electrograms
[uniDistInfo, uniDistVarnames, uniDistData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() 'Wave_uni_distal.csv']);

%get the roving unipolar proximal electrograms
[uniProxInfo, uniProxVarnames, uniProxData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() 'Wave_uni_proximal.csv']);


%% Parse electrogram data
% By loading the automaps (not needed for Precision but might be needed for
% EnSiteX)
% TODO

%% Calculate annotation times

% these are all in samples
refTick_adj   = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'Ref Tick'))); % _adj because these already reflect the user adjustments
rovTick_adj   = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'Rov Tick 1')));
leftCurtain     = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'left curtain (ms)'))) / 1000 * rovInfo.sampleFreq;
rightCurtain    = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'right curtain (ms)'))) / 1000 * rovInfo.sampleFreq;

startTime_s       = str2double(rovData(:,strcmpi(rovVarnames, 'startTime (abs)'))); % this comes from the roving wave file
refTime_s         = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'refTime (abs)'))); % this comes from the mapping file
adjTime_ms        = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'adjTime (ms)'))); % from the mapping file
lat_ms            = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'LAT'))); % from the mapping file
adjTime_samples   = adjTime_ms / 1000 * rovInfo.sampleFreq;

% annot method time based
refAnnot_time = (refTime_s - startTime_s) * rovInfo.sampleFreq;
latAnnot_time = (refTime_s + lat_ms/1000) * rovInfo.sampleFreq;

% annot method tick based
refAnnot_tick = refTick_adj;
latAnnot_tick = rovTick_adj;

% adjust time YES
%now we need to do nothing to the tick marks
refAnnot_tick_adj = refAnnot_tick;
latAnnot_tick_adj = latAnnot_tick;

%but we need to ADD the adjust time to the time-based times
refAnnot_time_adj = refAnnot_time + adjTime_samples;
latAnnot_time_adj = latAnnot_time + adjTime_samples;

% adjust time NO
%now we need to subtract the adjust time (converted into samples) to the tick marks
refAnnot_tick_noadj = refAnnot_tick - adjTime_samples;
latAnnot_tick_noadj = latAnnot_tick - adjTime_samples;

%but we do not need to do anything to the time based times
refAnnot_time_noadj = refAnnot_time;
latAnnot_time_noadj = latAnnot_time;

%TODO now we have annotations using all the methods we can check them ***

%Finally save the desired annotations
annotMethod = 'tickbased';
adjustTimes = 'no';
if strcmpi(annotMethod, 'timebased') && strcmpi(adjustTimes, 'yes')
    refAnnot = refAnnot_time_adj;
    latAnnot = latAnnot_time_adj;
end
if strcmpi(annotMethod, 'timebased') && strcmpi(adjustTimes, 'no')
    refAnnot = refAnnot_time_noadj;
    latAnnot = latAnnot_time_noadj;
end
if strcmpi(annotMethod, 'tickbased') && strcmpi(adjustTimes, 'yes')
    refAnnot = refAnnot_tick_adj;
    refAnnot = refAnnot_tick_adj;
end
if strcmpi(annotMethod, 'tickbased') && strcmpi(adjustTimes, 'no')
    refAnnot = refAnnot_tick_noadj;
    latAnnot = latAnnot_tick_noadj;
end

%% Save all files in the OpenEP structure

% General data
userdata = openep_createuserdata();
userdata.systemName = 'ensitex';
userdata.notes{1} = [date() ': Created'];
userdata.ensitexFolder = studyDir;
userdata.electric.sampleFrequency = rovInfo.sampleFreq;

% Geometry
userdata.surface.triRep = t;
surfaceData = data_geometry.dxgeo.surface_of_origin;
userdata = setSurfaceProperty(userdata, 'name', 'surfaceOfOrigin', 'map', surfaceData, 'definedOn', 'elements');
userdata.surface.normals = normals;

% Surface maps, removing invalid data beyond interpolation distance
userdata.surface.act_bip = [act bip];
userdata.surface.uni_imp_frc = [uni imp frc];

% Electric data
userdata.electric.electrodeNames_bip    = mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'Rov trace'));
userdata.electric.egmX                  = [str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'roving x'))) ...
                                           str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'roving y'))) ...
                                           str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'roving z')))];
userdata.electric.egmSurfX              = [str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'surface x'))) ...
                                           str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'surface y'))) ...
                                           str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'surface z')))];
userdata.electric.egmRef                = local_concatdata(refData(:,strcmpi(refVarnames,'signals')) ...
                                            ,refData(:,strcmpi(refVarnames,'Freeze Grp #')) ... 
                                            ,rovData(:,strcmpi(rovVarnames,'Freeze Grp #')));
userdata.electric.egm                   = local_concatdata(rovData(:,strcmpi(rovVarnames,'signals')),[],[]);

userdata.electric.annotations.referenceAnnot = refAnnot;
userdata.electric.annotations.mapAnnot  = latAnnot;
userdata.electric.annotations.woi       = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'left curtain (ms)'))) / 1000 * userdata.electric.sampleFrequency;
userdata.electric.annotations.woi(:,2)  = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'right curtain (ms)'))) / 1000 * userdata.electric.sampleFrequency;

userdata.electric.voltages.bipolar      = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'P-P')));
userdata.electric.include               = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'utilized')));
userdata.electric.names                 = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, '(Point #)')));

userdata.electric.electrodeNames_uni    = local_parseuninames(mappingPoints.data(:,strcmpi(mappingPoints.varnames,'Electrodes')));
userdata.electric.egmUniX               = userdata.electric.egmX; % note that we do not have co-ordinates for the second unipole
userdata.electric.egmUniSurfX           = userdata.electric.egmSurfX; % note that we do not have co-ordinates for the second unipole
userdata.electric.egmUni                = local_concatdata(uniDistData(:,strcmpi(uniDistVarnames,'signals')),[],[]);
userdata.electric.egmUni(:,:,2)         = local_concatdata(uniProxData(:,strcmpi(uniProxVarnames,'signals')),[],[]);

% we don't have the unipolar peak to peak voltages so we have to calculate them
userdata.electric.voltages.unipolar     = calculatePeak2PeakVoltage( userdata.electric.egmUni, userdata.electric.annotations.referenceAnnot, userdata.electric.annotations.woi );

% Temp - remove signalMaps which, if empty, prevents the file being loaded in EP Workbench
userdata.surface = rmfield(userdata.surface, 'signalMaps');
userdata.electric.tags = cell(length(userdata.electric.names),1);

% Encourage user to save the data
matFileFullPath = [];
if ~isempty(saveFileName_cli)
    save(saveFileName_cli, 'userdata');
    matFileFullPath = saveFileName_cli;
else
    defaultName = [mappingPoints.info.study '_' mappingPoints.info.mapName];
    defaultName(isspace(defaultName)) = '_';
    originalDir = cd();
    matFileFullPath = fullfile(saveDir, defaultName); %default
    cd(saveDir);
    [filename,saveDir] = uiputfile('*.mat', 'Save the userdata to disc for future rapid access?',defaultName);
    cd(originalDir);
    % We save as -v7 because it's faster to load in OpenEP-py than -v7.3,
    % and the saved file is significantly smaller compared to -v6 files.
    if filename ~= 0
        save([saveDir filename], 'userdata','-v7');
        matFileFullPath = fullfile(saveDir, filename);
    end
end

%% Local functions
    function pathName = local_findDirectory(stub, studyDir)
        allSubFolders = nameFolds(studyDir);
        thisFolder = allSubFolders(strstartcmpi(stub, allSubFolders));
        % Check if more than one folder meets the critiera, and ask the user
        % to choose
        if numel(thisFolder)>1
            warning(['IMPORTPRECISION_OPENEP: More than one candidate folder selected for the export of ***' stub '*** data. Please choose one folder ...'])
            [indx, tf] = listdlg('ListString', thisFolder ...
                ,'ListSize', [480 300] ...
                , 'name', ['Which is the correct ***' stub '*** folder?'] ...
                , 'selectionmode', 'single' ...
                );
            if ~tf
                error('IMPORTPRECISION_OPENEP: Operation cancelled')
            else
                thisFolder = thisFolder{indx};
            end
        end
        pathName = fullfile(studyDir, thisFolder);
    end

    function hd = local_homedirec()
        %HOMEDIREC returns the user's home directory.

        if ispc
            hd = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
        else
            hd = getenv('HOME');
        end
    end

    function matrixData = local_concatdata(cellData, freezeGroupIn, freezeGroupOut)
        % This function concatenates cell data into a matrix, optionally
        % based on the ordering specified by freezeGroupIn and freezeGroupOut
        f = waitbar(0, 'Reorganising data');
        if isempty(freezeGroupIn)
            repMatrix = 1:numel(cellData);
        else
            freezeGroupIn = str2double(freezeGroupIn);
            freezeGroupOut = str2double(freezeGroupOut);
            for iGrp = 1:numel(freezeGroupOut)
                repMatrix(iGrp,1) = find(freezeGroupIn==freezeGroupOut(iGrp));
            end
        end
        cellDataNew = cellData(repMatrix);
        nCell = numel(cellDataNew);
        matrixData = cellDataNew{1};
        for iCell = 2:nCell
            matrixData = [matrixData; cellDataNew{iCell}];
            waitbar(iCell/nCell, f);
        end

        % destroy the waitbar
        close(f)
    end

    function uniNames = local_parseuninames(A)
        % This function creates an Nx2 cell array for storing the unipole
        % names
        nPairs = size(A,1);
        uniNames = cell(nPairs,2);
        for iPair = 1:nPairs
            splt = strsplit(A{iPair});
            uniNames{iPair,1} = splt{1};
            uniNames{iPair,2} = splt{2};
        end
    end

end