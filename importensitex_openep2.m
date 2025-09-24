function [userdata, matFileFullPath] = importensitex_openep2(varargin)
% IMPORTENSITEX_OPENEP provides a data structure from multiple Precision files.
% Usage:
%   userdata = importprecision_openep(userinput)
%   userdata = importprecision_openep()
%   [userdata, matFileFullPath] = ...
% Where:
%   dirName is the directory with all of the files corresponding to a map
%   userdata is a single data structure
%   matFileFullPath is the path to the .mat file, if opened or saved
%
% IMPORTENSITEX_OPENEP accepts the following parameter-value pairs:
%   'savefilename'      {''}|string
%       The full path to the location in which to save the output.
%   'type'              {'standard'}|'omnipolar'
%       Specifies whether to import a standard map, based on local
%       activation time (LAT), bipolar (BIP) and unipolar (UNI) folders, or
%       whether to important an omnipolar dataset
%   'loadallwavefiles'  {false}|true
%       If set to true, the programme will load any wave files that are
%       extra to rov, ref and uni files into the ecg part of the userdata
%       structure (userdata.electric.ecg). This section of the data format
%       does not have its own positional information, unlike
%       userdata.electric.egm and userdata.electric.egmUni.
%
% IMPORTENSITEX_OPENEP is for parsing data from the EnsiteX mapping system.
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
%
% Author: Steven Williams (2025)
% SPDX-License-Identifier: Apache-2.0
%
% ---------------------------------------------------------------
% testing
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

%% Identify the folder containing the exported files
% Folder identification is either via the command line or via a pop up
% dialog box
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
% Additional command line input is parsed to determine the save location
% and whether a conventional or an omnipolar map is being assessed.
nStandardArgs = 1; % UPDATE VALUE
saveFileName_cli = '';
type = 'standard';
loadallwavefiles = 'false';
if nargin > nStandardArgs
    for i = nStandardArgs+1:2:nargin
        switch varargin{i}
            case 'savefilename'
                saveFileName_cli = varargin{i+1};
            case 'type'
                type = varargin{i+1};
            case 'loadallwavefiles'
                loadallwavefiles = varargin{i+1};
            otherwise
                error('IMPORTCARTO_MEM: Unrecognized input.')
        end
    end
end










%% Identify the relevant subfolders

% Naming of these folders needs to be done by the user either at the time
% of data export or afterwards. Although this is a manual step it avoids
% any ambiguity over which folder to import. Details are given in the SOP,
% "Instructions to convert Abbott Precision and EnSiteX data into OpenEP
% format"
switch type
    case 'standard'
        disp('IMPORTENSITEX_OPENEP: Type is standard');
        latMapDir = local_findDirectory('OpenEP_LAT', studyDir);
        ppMapDir = local_findDirectory('OpenEP_BIP', studyDir);
        uniMapDir = local_findDirectory('OpenEP_UNI', studyDir);
    case 'omnipolar'
        disp('IMPORTENSITEX_OPENEP: Type is omnipolar');
        latMapDir = local_findDirectory('OpenEP_OMNI_LAT', studyDir);
        ppMapDir = local_findDirectory('OpenEP_OMNI_PP', studyDir);
        uniMapDir = []; % TODO: is there a corollary of a uni map under omnipole polarity?
end

% Include the set of all the automap folders. Note that these are not yet
% being procesed
allSubFolds = nameFolds(studyDir);
automapDirs = allSubFolds(strstartcmpi('OpenEP_AutoMap', allSubFolds));










%% Parse the geometry and surface mapping data
% By loading the LAT XML file to get the geometry and the LAT map. In
% standard mode also load the PP and UNI maps to get substrate data. In
% omnipolar mode, we have the only LAT and PP maps (for now).
% TODO: check what is an UNI map recorded in omnipolar polarity.
data_geometry = loadprecision_modelgroups(fullfile(latMapDir{:}, 'Contact_Mapping_Model.xml'));
TRI = data_geometry.dxgeo.triangles;
X = data_geometry.dxgeo.vertices(:,1);
Y = data_geometry.dxgeo.vertices(:,2);
Z = data_geometry.dxgeo.vertices(:,3);
tr = TriRep(TRI, X, Y, Z);
t.X = tr.X;
t.Triangulation = tr.Triangulation;
normals = data_geometry.dxgeo.normals;

data_latMap = data_geometry; % no need to reload this
data_ppMap = loadprecision_modelgroups(fullfile(ppMapDir{:}, 'Contact_Mapping_Model.xml'));

act = data_latMap.dxgeo.act;
bip = data_ppMap.dxgeo.bip;

if ~isempty(uniMapDir)
    data_uniMap = loadprecision_modelgroups(fullfile(uniMapDir{:}, 'Contact_Mapping_Model.xml'));
    uni = data_uniMap.dxgeo.bip;
else
    uni = NaN(size(act));
end
imp = NaN(size(uni));
frc = NaN(size(uni));

%act(data_latMap.dxgeo.map_status==2) = NaN;
%bip(data_bipolarVoltageMap.dxgeo.map_status==2) = NaN;
%uni(data_unipolarVoltageMap.dxgeo.map_status==2) = NaN;










%% Parse electrogram data by loading the wave files (new equivalent of DxL files)

% work out the mapping points file names
switch type
    case 'standard'
        mapCSV = 'Map_LAT_bi.csv';
        voltCSV = 'Map_PP_bi.csv';

    case 'omnipolar'
        mapCSV = 'Map_LAT_omni.csv';
        voltCSV = 'Map_PP_omni.csv';

end
%load the mapping points data
[info, varnames, data] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() mapCSV]);
mappingPoints.info = info;
mappingPoints.varnames = varnames;
mappingPoints.data = data;

%additionally get the substrate mapping data for each point (there is unavoidable redundancy here)
[info, varnames, data] = loadensitex_dxldata([ppMapDir{:} filesep() 'Contact_Mapping' filesep() voltCSV]);
voltageData.info = info;
voltageData.varnames = varnames;
voltageData.data = data;

% do some simple checks for compatibility between the PP and LAT datasets
isError = false;
if mappingPoints.info.numPoints ~= voltageData.info.numPoints
    warning('IMPORTENSITEX_OPENEP: Mismatch between number of points in the voltage and activation time datasets');
    isError = true;
end
if ~strcmpi(mappingPoints.info.mapName, voltageData.info.mapName)
    warning('IMPORTENSITEX_OPENEP: Mismatch between map names in the voltage and activation time datasets');
    isError = true;
end
if ~strcmpi(mappingPoints.info.study, voltageData.info.study)
    warning('IMPORTENSITEX_OPENEP: Mismatch between study names in the voltage and activation time datasets');
    isError = true;
end
if isError
    error('IMPORTENSITEX_OPENEP: Error parsing data. See warnings above for hints');
end
% TODO: there are likely to be other checks we could add in here

% access the voltage data from the PP data and save along with the LAT data
ppStr = 'P-P'; ppValidStr = 'P-P valid';
mappingPoints.varnames{end+1} = ppStr;
mappingPoints.varnames{end+1} = ppValidStr;
ppData = voltageData.data(:,strcmpi(voltageData.varnames, ppStr));
ppValidData = voltageData.data(:,strcmpi(voltageData.varnames, ppValidStr));

% concatenate
mappingPoints.data = [mappingPoints.data ppData ppValidData];










%% Now load the electrogram data

% first load the rovinig trace (Wave_rov.csv),
% next load the reference trace (Wave_ref.csv),
% then load the unipolar electrograms, (Wave_uni_distal.csv, Wave_uni_along.csv), and
% finally load any other wave files that are present

%create loadedFiles boolean array to keep track of which wave files have already been loaded
allFiles = nameFiles([latMapDir{:} filesep() 'Contact_Mapping']);
waveFiles = find(~cellfun('isempty', regexp(allFiles, '^Wave', 'once')));
loadedFiles = false(size(waveFiles));

%get the reference electrograms
thisFilename = 'Wave_refs.csv';
[refInfo, refVarnames, refData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
loadedFiles(iThisFile) = true;


%get the roving bipolar electrograms
thisFilename = 'Wave_rov.csv';
[rovInfo, rovVarnames, rovData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
loadedFiles(iThisFile) = true;

%import the unipolar electrograms
switch type
    case 'standard'
        %import uni distal
        thisFilename = 'Wave_uni_distal.csv'
        [uniDistInfo, uniDistVarnames, uniDistData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
        iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
        loadedFiles(iThisFile) = true;

        %import uni proximal
        thisFilename = 'Wave_uni_proximal.csv';
        [uniProxInfo, uniProxVarnames, uniProxData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
        iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
        loadedFiles(iThisFile) = true;

    case 'omnipolar'
        %import uni across
        thisFilename = 'Wave_uni_across.csv';
        [uniAcrossInfo, uniAcrossVarnames, uniAcrossData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
        iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
        loadedFiles(iThisFile) = true;

        %import uni along
        thisFilename = 'Wave_uni_along.csv';
        [uniAlongInfo, uniAlongVarnames, uniAlongData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
        iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
        loadedFiles(iThisFile) = true;

        %import uni corner
        thisFilename = 'Wave_uni_corner.csv';
        [uniCornerInfo, uniCornerVarnames, uniCornerData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
        iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
        loadedFiles(iThisFile) = true;
end

%check for any other wave or map files
if any(~loadedFiles)
    % we have additional wave files, check whether to just load these or
    % ask the user what to do
    if loadallwavefiles
        extraFilesToLoad = allFiles(~loadedFiles);
        for iFile = 1:numel(extraFilesToLoad)
            if strcmpi(extraFilesToLoad{iFile}, 'Wave_refs2.csv')
                % This is to make sure that we do not attempt to load a
                % Wave_refs2.csv file, which for now seemt to be empty.
                % TODO: revisit Wave_refs2.csv files in the future if new
                % data is present
                continue
            else
                [extraFilesInfo{iFile}, extraFilesVarnames{iFile}, extraFilesData{iFile}] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() extraFilesToLoad{iFile}]);
            end
        end
    else
        warning('IMPORTENSITEX_OPENEP: Extra wave files identified, but code to ask the user what to do has not yet been implemented. For now if you want access to these wavefiles, re-run this programme with the option loadallwavefiles set to TRUE');
    end
end










%% Parse electrogram data
% By loading the automaps
% This may be necessary to give access to other data such as body surface
% ECG data
% TODO: Implement code to load and parse the automap data.










%% Calculate annotation times

% these are all in samples
refTick_adj   = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'Ref Tick'))); % _adj because these already reflect the user adjustments
rovTick_adj   = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'Rov Tick 1')));

%TODO: CHECK THESE TIMES ARE ADJUSTED APPROPRIATELY
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

userdata.electric.egmUni = zeros([size(userdata.electric.egm) size(userdata.electric.egmUniX,3)]);





% Unipolar electrogram locations are stored differently for standard and omnipolar configurations
switch type
    case 'standard'
        disp('IMPORTENSITEX_OPENEP: Parsing unipolar co-ordinates for omnipolar configuration ...');

        warning('IMPORTENSITEX_OPENEP: When a map is exported in non-unipolar mode we are not given the individual unipole co-ordinates ... assuming uni distal and uni proximal are at the same location')

        userdata.electric.egmUniX = [userdata.electric.egmX; userdata.electric.egmX];

        disp('IMPORTENSITEX_OPENEP: Finished parsing unipolar co-ordinates ...');
        
    case 'omnipolar'
        disp('IMPORTENSITEX_OPENEP: Parsing unipolar co-ordinates for omnipolar configuration ...');

        % now we need to locate the correct co-ordinates from mappingPoints,
        % without assuming that the order is correct. The order is however
        % usually (1) corner, (2) along, (3) across, so that is what we will
        % check first
        for iPoint = 1:size(userdata.electric.electrodeNames_uni,1)
            if strcmpi(userdata.electric.electrodeNames_uni{iPoint,1}, mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_Corner_Elec')))
                X = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_CornerX'));
                Y = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_CornerY'));
                Z = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_CornerZ'));
                userdata.electric.egmUniX(iPoint,1:3,1) = str2double([X Y Z]);

                % Store the relevant electrograms
                userdata.electric.egmUni(:,:,1) = local_concatdata(uniCornerData(:,strcmpi(uniCornerVarnames,'signals')),[],[]);
            else
                error('IMPORTENSITEX_OPENEP: Electrode naming mismatch')
            end
            if strcmpi(userdata.electric.electrodeNames_uni{iPoint,2}, mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_Along_Elec')))
                X = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_AlongX'));
                Y = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_AlongY'));
                Z = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_AlongZ'));
                userdata.electric.egmUniX(iPoint,1:3,2) = str2double([X Y Z]);

                % Store the relevant electrograms
                userdata.electric.egmUni(:,:,2) = local_concatdata(uniAlongData(:,strcmpi(uniAlongVarnames,'signals')),[],[]);
            else
                error('IMPORTENSITEX_OPENEP: Electrode naming mismatch')
            end
            if strcmpi(userdata.electric.electrodeNames_uni{iPoint,3}, mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_Across_Elec')))
                X = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_AcrossX'));
                Y = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_AcrossY'));
                Z = mappingPoints.data(iPoint,strcmpi(mappingPoints.varnames,'Uni_AcrossZ'));
                userdata.electric.egmUniX(iPoint,1:3,3) = str2double([X Y Z]);

                % Store the relevant elecctrograms
                userdata.electric.egmUni(:,:,3) = local_concatdata(uniAcrossData(:,strcmpi(uniAcrossVarnames,'signals')),[],[]);
            else
                error('IMPORTENSITEX_OPENEP: Electrode naming mismatch')
            end
        end
        userdata.electric.egmUniSurfX           = userdata.electric.egmUniX; % note that we do not have co-ordinates for the second unipole
        disp('IMPORTENSITEX_OPENEP: Finished parsing unipolar co-ordinates ...');
end

% Store the unipole electrograms into userdata.electric.egmUni matrix


% Store any additional signals in the ecg array
disp('dealing with ECG electrograms')
userdata.electric.ecgNames = {};
if ~isempty(extraFilesInfo)
    % for speed first work out the dimensions and pre-populate
    fWait = waitbar(0, 'Storing additional ECG names');
    for iEF = 1:numel(extraFilesInfo)
        if any(strcmpi(extraFilesVarnames{iEF}, 'signals'))
            % Then, this extra file contains signal data - store these in
            % the ECG array.

            % First check for unique electrode names in this file
            electrodes = unique(extraFilesData{iEF}(:,1));

            % Remove any non-ASCII characters, leading or trailing spaces
            % and duplicate rows (cElectrodes for 'clean electrodes')
            cElectrodes = unique(cellfun(@(s) strtrim(regexprep(s, '[^\x00-\x7F]', '')), electrodes, 'UniformOutput', false)); 

            % Add electrode names to the ecgNames cell array
            userdata.electric.ecgNames = union(userdata.electric.ecgNames, cElectrodes);
        else
            % Then, this extra file does not contain signal data. If the 
            % file has not already been imported (we do not yet have a
            % check for this) then it is likely to be an additional mapping
            % file. Do nothing for the time being.
        end
        waitbar(iEF/numel(extraFilesInfo),fWait);
    end
    close(fWait);

    % prepopulate for speed
    userdata.electric.ecg = zeros([size(userdata.electric.egm) size(userdata.electric.ecgNames,1)]);

    fWait = waitbar(0, 'Storing additional ECG data');
    for iEF = 1:numel(extraFilesInfo)
        % check if this is a signals file
        if any(strcmpi(extraFilesVarnames{iEF}, 'signals'))        
            % Next we iterate through every signal and work out where to
            % put it in the ECG array.

            for jSg = 1:size(extraFilesData{iEF},1)
                thisSig = extraFilesData{iEF}(jSg,strcmpi(extraFilesVarnames{iEF},'signals'));
                thisName = extraFilesData{iEF}(jSg,strcmpi(extraFilesVarnames{iEF},'Trace'));

                % clean the name
                thisName = strtrim(regexprep(thisName, '[^\x00-\x7F]', ''));

                userdata.electric.ecg(jSg,:,strcmpi(userdata.electric.ecgNames, thisName)) = thisSig{:};
            end
        else
            % Then, this extra file does not contain signal data. If the 
            % file has not already been imported (we do not yet have a
            % check for this) then it is likely to be an additional mapping
            % file. Do nothing for the time being.
        end
        waitbar(iEF/numel(extraFilesInfo),fWait);
    end
    close(fWait);

end







% we don't have the unipolar peak to peak voltages so we have to calculate them
% userdata.electric.voltages.unipolar     = calculatePeak2PeakVoltage( userdata.electric.egmUni, userdata.electric.annotations.referenceAnnot, userdata.electric.annotations.woi );

% Temp - remove signalMaps which, if empty, prevents the file being loaded in EP Workbench
userdata.surface = rmfield(userdata.surface, 'signalMaps');
userdata.electric.tags = cell(length(userdata.electric.names),1);










%% Encourage user to save the data
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
        save([saveDir filename], 'userdata','-v7.3');
        matFileFullPath = fullfile(saveDir, filename);
    end
end










%% Local functions

    function [info, varnames, data] = local_concatenatedxldatasets(inInfo, inVarnames, inData)
        % LOCAL_CONCATENATEDXLDATASETS is a function which combines
        % datasets together into a single set of structures
        % How many datasets are there

        disp('IMPORTENSITEX_OPENEP: concatenating data using local_concatenatedxldatasets')
        numDatasets = numel(inInfo);

        % Concatenate the info data.
        info = inInfo{1};

        header = info.header;
        info.header = [];
        info.header{1} = header;

        fname = info.filename;
        info.filename = [];
        info.filename{1} = fname;

        for iD = 2:numDatasets % D for dataset
            info.numPoints = info.numPoints + inInfo{iD}.numPoints;
            info.header{iD} = inInfo{iD}.header;
            info.filename{iD} = inInfo{iD}.filename;
        end

        info.header = info.header';
        info.filename = info.filename';

        % TODO: Add checks to make sure no other important aspects of rovInfoA and rovInfoB change.

        % Concatenate the varnames data.
        % Save the first varnames, but check with the subsequent
        % varnames for any differences. Throw an error if found.
        varnames = inVarnames{1};
        for iD = 2:numDatasets
            if ~comparestructure(inVarnames{iD}, inVarnames{iD-1})
                error('IMPORTENSITEX_OPENEP: Problem with source data - roving variable names changes between wave_across and wave_along');
            end
        end

        % Concatenate the data
        data = inData{1};
        for iD = 2:numDatasets
            data = [data; inData{iD}];
        end
    end

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
        %uniNames = cell(nPairs,2);
        for iPair = 1:nPairs
            splt = strsplit(A{iPair});
            for j = 1:numel(splt)
                uniNames{iPair,j} = splt{j};
            end
            % uniNames{iPair,1} = splt{1};
            % uniNames{iPair,2} = splt{2};
        end
    end

end