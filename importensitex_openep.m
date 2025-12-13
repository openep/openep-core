function [userdata, matFileFullPath] = importensitex_openep(varargin)
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
%   'maptoread'         {''}|string|double
%       Specifies which map to read. Can be a string referring
%       to the map name or a double referring to the number of points in the
%       map. If there are multiple maps with the same number of points an error
%       will be thrown.
%       Specifies whether to import the bipolar, omnipolar or unipolar
%       electrograms. 
%   'egmtype'           'bi'|'omni'|'uni'
%   'maptype'           {'asegm'}|'all'
%       Specifies whether to only import the surface map linked to the
%       chosen electrograms; or to search for additional matching
%       Contact_Mapping_Model.xml files and import them too.
%   'savefilename'      {''}|string
%       The full path to the location in which to save the output.

%
% IMPORTENSITEX_OPENEP is for parsing data from the EnsiteX mapping system.
% The function handles data in bipolar, omnipolar and unipolar format. One,
% two or all of these formats can be present. The function works on a
% single map.
%
% Some important considerations:
%
% re: egmtype There is no option to import all the electrograms. - if this
% is desired multiple imports must be run creating different OpenEP files.
%
% (1) If maptype is 'all', then the following convetions apply:
%       - act and bip taken from bipolar map folder
%       - uni taken from unipolar map folder, imp and frc are not populated
%       - all other available maps are stored as surface properties
%       - if any of these folders are missing, a warning is given and the
%       relevant data fields are empty
%       - geometry is taken from the first available folder in the order of
%       preference of bipolar > omnipolar > unipolar
% (2) If maptype is one of 'bipolar', 'omnipolar' or 'unipolar' then not
%   all data fields will be populated. Specifically:
%       - 'bipolar'   - act, bip populated; uni not available
%       - 'omnipolar' - act populated, bip and uni not available
%       - 'unipolar'  - act and uni populated, bip not available
%       - all other available maps are stored as surface properties
%       - geometry is taken from the specified folder
%       - if the specified folder does not exist, an error is thrown
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
mapToRead = '';
egmtype = '';
maptype = 'asegm';
saveFileName = '';

if nargin > nStandardArgs
    for i = nStandardArgs+1:2:nargin
        switch varargin{i}
            case 'maptoread'
                mapToRead = varargin{i+1};
            case 'egmtype'
                egmtype = varargin{i+1};
            case 'maptype'
                maptype = varargin{i+1};
            case 'savefilename'
                saveFileName = varargin{i+1};
            otherwise
                error('IMPORTENSITEX_OPENEP: Unrecognized input.')
        end
    end
end










%% Identify all available maps and export styles (uni, omni, bip)

% This logic is reasonably robust but there are some requirements. 
% All wave and Map files that are related to each other are stored in 
% separate directories. This is the default way that the files come out of the 
% system, named by a timestamp plus any other text that the user added. 
% However, if the user has moved these files 
% to a different location this code will throw an error since _almost 
% certainly_ the number of points in wave and map files will no longer 
% match. 

csvFiles = local_findAllCsvFiles(studyDir);
allCsvHeaders = [];
for iCsv = 1:numel(csvFiles)
    csvHeader = local_loadCsvFileHeader(csvFiles{iCsv});
    if ~isempty(csvHeader)
        allCsvHeaders{end+1} = csvHeader;
    end
end

% Identify the names of unique maps in the dataset
uniqueMapNames = unique( cellfun(@(s) s.mapName, allCsvHeaders, 'UniformOutput', false) );

% Identify all folder locations in which files pertaining to each map are stored
numMaps = numel(uniqueMapNames);
locations = cell(numel(uniqueMapNames), 2);
for iMap = 1:numMaps
    name = uniqueMapNames{iMap};

    % Select structures belonging to this map name
    idx = cellfun(@(s) strcmp(s.mapName, name), allCsvHeaders);

    % Extract filenames for this map
    files = cellfun(@(s) s.filename, allCsvHeaders(idx), 'UniformOutput', false);

    % Extract unique folder paths
    folders = unique(cellfun(@fileparts, files, 'UniformOutput', false));

    % Store results
    locations{iMap,1} = strrep(name, sprintf('\t'), ' ');
    locations{iMap,2} = folders;
end

% CHECK1A: Check that all map files within each directory pertaining to an 
% individual map have the same type (uni, omni, bi). There is redundancy in
% this section since we re-read the CSV file headers per folder, and we
% have already read them all collectively. However, this avoids any
% ambiguity about what files are being read and compared.
% CHECK2: Check that all wave files within each directory pertaining to an
% individual map have the same number of electrograms as number of points
% in the map files; with the one exception being wave_rov files
allMapTypes = {};
for iMap = 1:numMaps
    folders = locations{iMap,2};
    numFolds = numel(folders);

    folderMapType = {};
    numberOfPoints = [];
    for iFolder = 1:numFolds

        % identify the CSV files for this map
        csvFilesInThisFolder = local_findAllCsvFiles(locations{iMap,2}{iFolder});

        % read the header of these CSV files
        csvFilesInThisFolderHeaders = [];
        for iCsv = 1:numel(csvFilesInThisFolder)
            csvHeader = local_loadCsvFileHeader(csvFilesInThisFolder{iCsv});
            if ~isempty(csvHeader)
                csvFilesInThisFolderHeaders{end+1} = csvHeader;
            end
        end

        % get all the map types for these csv files
        allMapTypesInThisFolder = {};
        mappingCsvFiles = [];
        for iCsv = 1:numel(csvFilesInThisFolderHeaders)
            thisMapType = csvFilesInThisFolderHeaders{iCsv}.mapType;
            thisMapNumPts = csvFilesInThisFolderHeaders{iCsv}.numPoints;
            if ~strcmp(thisMapType, 'N/A')
                allMapTypesInThisFolder{end+1} = thisMapType;
                mappingCsvFiles{end+1} = csvFilesInThisFolderHeaders{iCsv}; % this ensure we are only dealing with mapping files, and not wave files
            end
        end

        [mapTypesAreTheSame, identifiedType] = local_areMapSuffixesUniform(allMapTypesInThisFolder);

        [numberOfMappingPointsAreTheSame, numMappingPts] = local_areMapPointNumbersUniform(mappingCsvFiles);

        % check that the map types are the same
        if mapTypesAreTheSame
            disp(['CHECK1 TEST PASSED: map types check passed for map ' num2str(iMap) ' folder ' num2str(iFolder)]);
            
            % store the type of maps that are in this folder
            folderMapType{end+1} = identifiedType;
        else
            error(['CHECK1 TEST FAILED: map types check failed for map ' num2str(iMap) ' folder ' num2str(iFolder)])
        end

        % check that the number of mapping points are the same
        if numberOfMappingPointsAreTheSame
            disp(['CHECK2 TEST PASSED: number of mapping points are the same for map ' num2str(iMap) ' folder ' num2str(iFolder)]);

            % store the number of mapping points
            numberOfPoints(end+1) = numMappingPts;
        else
            error(['CHECK2 TEST FAILED: number of mapping points are different for map ' num2str(iMap) ' folder ' num2str(iFolder)]);
        end
    end
    locations{iMap,3} = folderMapType;
    locations{iMap,5} = numberOfPoints;
end

% Identify the geoemtry location (Contact_Mapping_Model.xml)
% Option A - this file is 1 level up from the mapping/wave files
% Option B - this file is within the same folder as the mapping/wave files
% Option C - give an error
numMaps = numel(locations(:,1));
for iMap = 1:numMaps
    numFolders = numel(locations{1,2});
    for jFolder = 1:numFolders
        thisFolderPath = locations{iMap, 2}{jFolder};
        [path, ~] = fileparts(thisFolderPath);
        geomFile = [path filesep() 'Contact_Mapping_Model.xml'];
        if isfile(geomFile)
            locations{iMap,4}{jFolder} = geomFile;
        else
            geomFile = [thisFolderPath filesep() 'Contact_Mapping_Model.xml'];
            if isfile(geomFile)
                locations{iMap,4}{jFolder} = geomFile;
            else
                error('IMPORTENSITEX_OPENEP: Unable to idetnify the Contact_Mapping_Model.xml file');
            end
        end
    end
end






%% Save all the data we have worked out in a table for easy access
variableNames = {'mapname', 'egmfiles', 'egmtype', 'mapfiles', 'numpts'};
T = cell2table(locations, 'variablenames', variableNames);








%% Identify the relevant subfolders

% Naming of these folders needs to be done by the user either at the time
% of data export or afterwards. Although this is a manual step it avoids
% any ambiguity over which folder to import. Details are given in the SOP,
% "Instructions to convert Abbott Precision and EnSiteX data into OpenEP
% format"
% omniDir = local_findDirectory('omnipole', studyDir);
% bipDir = local_findDirectory('bipole', studyDir);
% uniDir = local_findDirectory('unipole', studyDir);






%% Ask the user which map they want to import
names = T.mapname;
numPtsPerMap = T.numpts;

if isempty(mapToRead)
    [selection,ok] = listdlg(     'ListString', names ...
        , 'SelectionMode', 'single' ...
        , 'PromptString', 'Which map do you want to access?' ...
        , 'ListSize', [300 300] ...
        );
    if ~ok
        return
    end
    mapToRead = names{selection};
else
    if isnumeric(mapToRead)
        selection = numel(find(numPtsPerMap==mapToRead));
        if numel(selection)>1
            error(['IMPORTENSITEX_OPENEP: Multiple maps with ' ...
                num2str(mapToRead) ...
                ' points identified. Use an alternative method to identify map.']);
        elseif isempty(selection)
            error(['IMPORTENSITEX_OPENEP: No map with ' ...
                num2str(mapToRead) ...
                ' points identified. Check the number of points specified is correct.']);
        else
            [selection, ~] = find(numPtsPerMap==mapToRead);
        end
    elseif ischar(mapToRead)
        selection = find(strstartcmpi(mapToRead, names));
    end
end
mapID = selection; % calling it map ID to be more understandable. MapID maps into rows of T.







%% Ask the user which mapping style they want to import, based on the available mapping styles
names = T(mapID, 'egmtype');
names = names{1,:};
uNames = unique(names); % convert to cell array and identify the unique names
if isempty(egmtype)
    [selection,ok] = listdlg(     'ListString', uNames ...
        , 'SelectionMode', 'single' ...
        , 'PromptString', 'Which electrogram type do you want?' ...
        , 'ListSize', [300 300] ...
        );
    if ~ok
        return
    end
end
reqEgmType = names{selection};
egmID = find(strcmpi(names, reqEgmType)); 
% note that egmID by itself is not interpretable, but it indexes into T
% table entries to ensure that the desired electrograms are read
egmtype = reqEgmType;

if numel(egmID)>1
    warningMessage = ['Multiple ' reqEgmType ' electrograms identified for map ' mapToRead '. Which folder of electrograms do you want to import?'];
    warning(['IMPORTENSITEX_OPENEP: ' warningMessage]);

    names = T.egmfiles(1,egmID);
    shortNames = local_lastTwoParts(names);

    [selection,ok] = listdlg(     'ListString', shortNames ...
        , 'SelectionMode', 'single' ...
        , 'PromptString', warningMessage ...
        , 'ListSize', [600 300] ...
        );
    if ~ok
        return
    end
    egmID = egmID(selection);
    egmtype = shortNames{selection};
end


    





%% Parse the geometry and surface mapping data
% By loading the relevant Contact_Mapping_Model XML file to get the geometry

contactMappingModel = T(mapID,:).mapfiles{egmID};
data_geometry = loadprecision_modelgroups(contactMappingModel);

% switch maptype
%     case 'bipolar'
%         data_geometry = loadprecision_modelgroups(fullfile(bipDir, 'Contact_Mapping_Model.xml'));
% 
%     case 'omnipolar'
%         data_geometry = loadprecision_modelgroups(fullfile(omniDir, 'Contact_Mapping_Model.xml'));
% 
%     case 'unipolar'
%         data_geometry = loadprecision_modelgroups(fullfile(uniDir, 'Contact_Mapping_Model.xml'));
% 
%     case 'all'
%         if isfolder(bipDir)
%             data_geometry = loadprecision_modelgroups(fullfile(bipDir, 'Contact_Mapping_Model.xml'));
%         elseif isfolder(omniDir)
%             data_geometry = loadprecision_modelgroups(fullfile(omniDir, 'Contact_Mapping_Model.xml'));
%         elseif isfolder(uniDir)
%             data_geometry = loadprecision_modelgroups(fullfile(uniDir, 'Contact_Mapping_Model.xml'));
%         end
% end

TRI = data_geometry.dxgeo.triangles;
X = data_geometry.dxgeo.vertices(:,1);
Y = data_geometry.dxgeo.vertices(:,2);
Z = data_geometry.dxgeo.vertices(:,3);
tr = TriRep(TRI, X, Y, Z);
t.X = tr.X;
t.Triangulation = tr.Triangulation;
normals = data_geometry.dxgeo.normals;






%% Parse the mapping data according to the users wishes
% Note that in this section, any time we store mapping data we also must
% check the map status to determine whether values should be replaced by
% NaN values.
lenX = size(tr.X,1);
act = NaN(lenX,1);
bip = NaN(lenX,1);
uni = NaN(lenX,1);
mapData = [];
mapType = [];

actSaved = false;
bipSaved = false;
uniSaved = false;

switch maptype
    case 'asegm'
        if ~isempty(data_geometry.dxgeo.act)
            act = data_geometry.dxgeo.act;
            iStatus = data_geometry.dxgeo.map_status;
            act(iStatus==2) = NaN;
            actSaved = true;
        end
        if ~isempty(data_geometry.dxgeo.bip)
            bip = data_geometry.dxgeo.bip;
            iStatus = data_geometry.dxgeo.map_status;
            bip(iStatus==2) = NaN;
            bipSaved = true;
        end
        if ~isempty(data_geometry.dxgeo.uni)
            uni = data_geometry.dxgeo.uni;
            iStatus = data_geometry.dxgeo.map_status;
            uni(iStatus==2) = NaN;
            uniSaved = true;
        end
        if isfield(data_geometry.dxgeo, 'mapdata')
            if ~isempty(data_geometry.dxgeo.mapdata)
                mapData = data_geometry.dxgeo.mapdata;
                mapType = data_geometry.dxgeo.maptype;
                iStatus = data_geometry.dxgeo.map_status;
                mapData(iStatus==2) = NaN;
            end
        end

    case 'all'

        % Lots of logic has to go into here - finding all XML files in
        % folders or subfolders, loading these XML files, checking whether
        % the geometry matches, if it does, load the corresponding map into
        % the right place (act, bip, uni or mapData), removing values that
        % should be NaN along the way.

        % First find all XML files in folder or subfolders
        xmlFiles = local_findAllXmlFiles(studyDir);

        % Load all these XML files
        for iXml = 1:numel(xmlFiles)
            dataXml{iXml} = loadprecision_modelgroups(xmlFiles{iXml});
        end

        % Compare the geometry between the XML files and the existing geometry
        % We define a match as an exact match of vetcies, triangles and
        % normals.
        for iXml = 1:numel(xmlFiles)
            fileIsValid(iXml) = local_compareXmlFiles(dataXml{iXml}, data_geometry);
        end

        % For every XML file that has a matching geometry, load the corresponding map
        for iXml = 1:numel(xmlFiles)
            dataIdentified = false;
            if fileIsValid(iXml)
                % First check for any of act, bip or uni
                if ~isempty(dataXml{iXml}.dxgeo.act)
                    if ~actSaved
                        act = dataXml{iXml}.dxgeo.act;
                        iStatus = dataXml{iXml}.dxgeo.map_status;
                        act(iStatus==2) = NaN;
                        actSaved = true;
                    else
                        warning(['IMPORTENSITEX_OPENEP: Multiple local activation time surface maps identified. ...' ...
                            'The first identified map comes from the file ', dataXml{iXml}.fileLoaded, ...
                            ' and is stored in .act_bip. The remaining maps are stored as surface properties.']);
                        mapData{end+1} = dataXml{iXml}.dxgeo.act;
                        mapType{end+1} = ['Additional LAT map ' num2str(numel(mapType))];
                        iStatus = dataXml{iXml}.dxgeo.map_status;
                        mapData{end}(iStatus==2) = NaN;
                    end
                    dataIdentified = true;

                end
                if ~isempty(dataXml{iXml}.dxgeo.bip)
                    if ~bipSaved
                        bip = dataXml{iXml}.bip;
                        iStatus = dataXml{iXml}.dxgeo.map_status;
                        bip(iStatus==2) = NaN;
                        bipSaved = true;
                    else
                        warning(['IMPORTENSITEX_OPENEP: Multiple bipolar voltage maps identified. ...' ...
                            'The first identified map comes from the file ', dataXml{iXml}.fileLoaded, ...
                            ' and is stored in .act_bip. The remaining maps are stored as surface properties.']);
                        mapData{end+1} = dataXml{iXml}.dxgeo.bip;
                        mapType{end+1} = ['Additional BIP map ' num2str(nunmel(mapType))];
                        iStatus = dataXml{iXml}.dxgeo.map_status;
                        mapData{end}(iStatus==2) = NaN;
                    end
                    dataIdentified = true;

                end
                if ~isempty(dataXml{iXml}.dxgeo.uni)
                    if ~uniSaved
                        uni = dataXml{iXml}.uni;
                        iStatus = dataXml{iXml}.dxgeo.map_status;
                        uni(iStatus==2) = NaN;
                        uniSaved = true;

                    else
                        warning(['IMPORTENSITEX_OPENEP: Multiple unipolar voltage maps identified. ...' ...
                            'The first identified map comes from the file ', dataXml{iXml}.fileLoaded, ...
                            ' and is stored in .uni_imp_frc. The remaining maps are stored as surface properties.']);
                        mapData{end+1} = dataXml{iXml}.dxgeo.uni;
                        mapType{end+1} = ['Additional UNI map ' num2str(nunmel(mapType))];
                        iStatus = dataXml{iXml}.dxgeo.map_status;
                        mapData{end}(iStatus==2) = NaN;
                    end
                    dataIdentified = true;

                end

                % Then check for any other mapping files
                if ~dataIdentified
                    mapData{end+1} = dataXml{iXml}.dxgeo.mapdata;
                    mapType{end+1} = dataXml{iXml}.dxgeo.maptype;

                    iStatus = dataXml{iXml}.dxgeo.map_status;
                    mapData{end}(iStatus==2) = NaN;

                end
            else
                warning(['IMPORTENSITEX_OPENEP: An XML mapping file which does ...' ...
                    'not match the loaded geometry has been identified. File ...' ...
                    , dataXml{iXml}.fileLoaded ' will be ignored.'])
                continue;

            end
        end
end

% switch maptype
    % case 'bipolar'
    %     % we know we have a bipolar map of some sort, so we will check for
    %     % an activation map, a voltage map or any other maps. We know we
    %     % will not have a unipolar map so we will set uni to [];
    %     act = data_geometry.dxgeo.act;
    %     bip = data_geometry.dxgeo.bip;
    %     uni = [];
    %     mapData = data_geometry.dxgeo.mapdata;
    %     mapType = data_geometry.dxgeo.maptype;
    % 
    %     iStatus = data_geometry.dxgeo.map_status;
    %     act(iStatus==2) = NaN;
    %     bip(iStatus==2) = NaN;
    %     mapData(iStatus==2) = NaN;
    % 
    % case 'omnipolar'
    %     % we know we will have an omnipolar map of some sort, but we will
    %     % not have a conventional bipolar LAT map, bipolar voltage map or
    %     % unipolar voltage map, so we will set act, uni and bip to [];
    %     act = [];
    %     bip = [];
    %     uni = [];
    %     mapData = data_geometry.dxgeo.mapdata;
    %     mapType = data_geometry.dxgeo.maptype;
    % 
    %     iStatus = data_geometry.dxgeo.map_status;
    %     mapData(iStatus==2) = NaN;
    % 
    % case 'unipolar'
    %     % we know we will have a unipolar map of some sort, but we will not
    %     % have a convetional bipolar LAT map, or bipolar votlage map, so we
    %     % will check for a uni voltage map and set act and bip to[];
    %     act = [];
    %     bip = [];
    %     uni = data_geometry.dxgeo.uni;
    %     mapData = data_geometry.dxgeo.mapdata;
    %     mapType = data_geometry.dxgeo.maptype;
    % 
    %     iStatus = data_geometry.dxgeo.map_status;
    %     uni(iStatus==2) = NaN;
    %     mapData(iStatus==2) = NaN;
% 
%     case 'all'
%         act = [];
%         bip = [];
%         uni = [];
%         mapData = [];
%         mapType = [];
% 
%         % Lots of logic has to go into here - finding all XML files in
%         % folders or subfolders, loading these XML files, checking whether
%         % the geometry matches, if it does, load the corresponding map into
%         % the right place (act, bip, uni or mapData), removing values that
%         % should be NaN along the way.
% 
%         % First find all XML files in folder or subfolders
%         xmlFiles = local_findAllXmlFiles(studyDir);
% 
%         % Load all these XML files
%         for iXml = 1:numel(xmlFiles)
%             dataXml{iXml} = loadprecision_modelgroups(xmlFiles{iXml});
%         end
% 
%         % Compare the geometry between the XML files and the existing geometry
%         % We define a match as an exact match of vetcies, triangles and
%         % normals.
%         for iXml = 1:numel(xmlFiles)
%             fileIsValid(iXml) = local_compareXmlFiles(dataXml{iXml}, data_geometry);
%         end
% 
%         % For every XML file that has a matching geometry, load the corresponding map
%         for iXml = 1:numel(xmlFiles)
%             dataIdentified = false;
%             if fileIsValid(iXml)
%                 % First check for any of act, bip or uni
%                 if ~isempty(dataXml{iXml}.dxgeo.act)
%                     if isempty(act)
%                         act = dataXml{iXml}.dxgeo.act;
% 
%                         iStatus = dataXml{iXml}.dxgeo.map_status;
%                         act(iStatus==2) = NaN;
% 
%                     else
%                         warning(['IMPORTENSITEX_OPENEP: Multiple local activation time surface maps identified. ...' ...
%                             'The first identified map comes from the file ', dataXml{iXml}.fileLoaded, ...
%                             ' and is stored in .act_bip. The remaining maps are stored as surface properties.']);
%                         mapData{end+1} = dataXml{iXml}.dxgeo.act;
%                         mapType{end+1} = ['Additional LAT map ' num2str(nunmel(mapType))];
% 
%                         iStatus = dataXml{iXml}.dxgeo.map_status;
%                         mapData{end}(iStatus==2) = NaN;
% 
%                     end
%                     dataIdentified = true;
% 
%                 end
%                 if ~isempty(dataXml{iXml}.dxgeo.bip)
%                     if isempty(bip)
%                         bip = dataXml{iXml}.bip;
% 
%                         iStatus = dataXml{iXml}.dxgeo.map_status;
%                         bip(iStatus==2) = NaN;
% 
%                     else
%                         warning(['IMPORTENSITEX_OPENEP: Multiple bipolar voltage maps identified. ...' ...
%                             'The first identified map comes from the file ', dataXml{iXml}.fileLoaded, ...
%                             ' and is stored in .act_bip. The remaining maps are stored as surface properties.']);
%                         mapData{end+1} = dataXml{iXml}.dxgeo.bip;
%                         mapType{end+1} = ['Additional BIP map ' num2str(nunmel(mapType))];
% 
%                         iStatus = dataXml{iXml}.dxgeo.map_status;
%                         mapData{end}(iStatus==2) = NaN;
% 
%                     end
%                     dataIdentified = true;
% 
%                 end
%                 if ~isempty(dataXml{iXml}.dxgeo.uni)
%                     if isempty(uni)
%                         uni = dataXml{iXml}.uni;
% 
%                         iStatus = dataXml{iXml}.dxgeo.map_status;
%                         uni(iStatus==2) = NaN;
% 
%                     else
%                         warning(['IMPORTENSITEX_OPENEP: Multiple unipolar voltage maps identified. ...' ...
%                             'The first identified map comes from the file ', dataXml{iXml}.fileLoaded, ...
%                             ' and is stored in .uni_imp_frc. The remaining maps are stored as surface properties.']);
%                         mapData{end+1} = dataXml{iXml}.dxgeo.uni;
%                         mapType{end+1} = ['Additional UNI map ' num2str(nunmel(mapType))];
% 
%                         iStatus = dataXml{iXml}.dxgeo.map_status;
%                         mapData{end}(iStatus==2) = NaN;
% 
%                     end
%                     dataIdentified = true;
% 
%                 end
% 
%                 % Then check for any other mapping files
%                 if ~dataIdentified
%                     mapData{end+1} = dataXml{iXml}.dxgeo.mapdata;
%                     mapType{end+1} = dataXml{iXml}.dxgeo.maptype;
% 
%                     iStatus = dataXml{iXml}.dxgeo.map_status;
%                     mapData{end}(iStatus==2) = NaN;
% 
%                 end
%             else
%                 warning(['IMPORTENSITEX_OPENEP: An XML mapping file which does ...' ...
%                     'not match the loaded geometry has been identified. File ...' ...
%                     , dataXml{iXml}.fileLoaded ' will be ignored.'])
%                 continue;
% 
%             end
%         end
% end

% IMP and FRC are not currently available through the EnsiteX export
% options

imp = NaN(size(uni));
frc = NaN(size(uni));

disp('!!!! FINISHED PARSING MAPPING DATA ACCCORDING TO USER WISHES !!!!')








%% Parse annotation metrics by loading the Map files

mappingPointsFolder = T(mapID,:).egmfiles{egmID};
csvFiles = local_findAllCsvFiles(mappingPointsFolder);
csvHeaders = [];
for iFile = 1:numel(csvFiles)
    tempCsvHeader = local_loadCsvFileHeader(csvFiles{iFile});
    if ~strcmp(tempCsvHeader.mapType, 'N/A') % we want the files where mapType is NOT N/A
        csvHeaders{end+1} = tempCsvHeader;
    end
end
for iFile = 1:numel(csvHeaders)
    [info, varnames, data] = loadensitex_dxldata(csvHeaders{iFile}.filename);
    mappingData{iFile}.info = info;
    mappingData{iFile}.varnames = varnames;
    mappingData{iFile}.data = data;
end


% work out the mapping points file names
% switch type
%     case 'standard'
%         mapCSV = 'Map_LAT_bi.csv';
%         voltCSV = 'Map_PP_bi.csv';
% 
%     case 'omnipolar'
%         mapCSV = 'Map_LAT_omni.csv';
%         voltCSV = 'Map_PP_omni.csv';
% 
% end
% %load the mapping points data
% latMapDir = []; %TEMP
% [info, varnames, data] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() mapCSV]);
% mappingPoints.info = info;
% mappingPoints.varnames = varnames;
% mappingPoints.data = data;
% 
% ppMapDir = []; % TEMP
% %additionally get the substrate mapping data for each point (there is unavoidable redundancy here)
% [info, varnames, data] = loadensitex_dxldata([ppMapDir{:} filesep() 'Contact_Mapping' filesep() voltCSV]);
% voltageData.info = info;
% voltageData.varnames = varnames;
% voltageData.data = data;
% 
% % do some simple checks for compatibility between the PP and LAT datasets
% isError = false;
% if mappingPoints.info.numPoints ~= voltageData.info.numPoints
%     warning('IMPORTENSITEX_OPENEP: Mismatch between number of points in the voltage and activation time datasets');
%     isError = true;
% end
% if ~strcmpi(mappingPoints.info.mapName, voltageData.info.mapName)
%     warning('IMPORTENSITEX_OPENEP: Mismatch between map names in the voltage and activation time datasets');
%     isError = true;
% end
% if ~strcmpi(mappingPoints.info.study, voltageData.info.study)
%     warning('IMPORTENSITEX_OPENEP: Mismatch between study names in the voltage and activation time datasets');
%     isError = true;
% end
% if isError
%     error('IMPORTENSITEX_OPENEP: Error parsing data. See warnings above for hints');
% end
% % TODO: there are likely to be other checks we could add in here
% 
% % access the voltage data from the PP data and save along with the LAT data
% ppStr = 'P-P'; ppValidStr = 'P-P valid';
% mappingPoints.varnames{end+1} = ppStr;
% mappingPoints.varnames{end+1} = ppValidStr;
% ppData = voltageData.data(:,strcmpi(voltageData.varnames, ppStr));
% ppValidData = voltageData.data(:,strcmpi(voltageData.varnames, ppValidStr));
% 
% % concatenate
% mappingPoints.data = [mappingPoints.data ppData ppValidData];





%% Parse electrogram data by loading the Wave files

wavesFolder = T(mapID,:).egmfiles{egmID};

%get the reference electrograms
thisFilename = 'Wave_refs.csv';
[refInfo, refVarnames, refData] = loadensitex_dxldata([wavesFolder filesep() thisFilename]);

%get the roving electrograms
thisFilename = 'Wave_rov.csv';
[rovInfo, rovVarnames, rovData] = loadensitex_dxldata([wavesFolder filesep() thisFilename]);

%import the unipolar electrograms
switch egmtype
    case 'bi'
        %import uni distal
        thisFilename = 'Wave_uni_distal.csv';
        [uniDistInfo, uniDistVarnames, uniDistData] = loadensitex_dxldata([wavesFolder filesep() thisFilename]);

        %import uni proximal
        thisFilename = 'Wave_uni_proximal.csv';
        [uniProxInfo, uniProxVarnames, uniProxData] = loadensitex_dxldata([wavesFolder filesep() thisFilename]);

    case 'omni'
        %import uni across
        thisFilename = 'Wave_uni_across.csv';
        [uniAcrossInfo, uniAcrossVarnames, uniAcrossData] = loadensitex_dxldata([wavesFolder filesep() thisFilename]);

        %import uni along
        thisFilename = 'Wave_uni_along.csv';
        [uniAlongInfo, uniAlongVarnames, uniAlongData] = loadensitex_dxldata([wavesFolder filesep() thisFilename]);

        %import uni corner
        thisFilename = 'Wave_uni_corner.csv';
        [uniCornerInfo, uniCornerVarnames, uniCornerData] = loadensitex_dxldata([wavesFolder filesep() thisFilename]);

    case 'uni'
        % In the case of uni electrogram mode, the roving electrogram is
        % the unipolar electrogram used to create the map and there are no
        % additional unipolar wave files so there is nothing else to do
end





% %% Now load the electrogram data by loading the Wave files (new equivalent of DxL files)
% 
% % first load the rovinig trace (Wave_rov.csv),
% % next load the reference trace (Wave_ref.csv),
% % then load the unipolar electrograms, (Wave_uni_distal.csv, Wave_uni_along.csv), and
% % finally load any other wave files that are present
% 
% %create loadedFiles boolean array to keep track of which wave files have already been loaded
% allFiles = nameFiles([latMapDir{:} filesep() 'Contact_Mapping']);
% waveFiles = find(~cellfun('isempty', regexp(allFiles, '^Wave', 'once')));
% loadedFiles = false(size(waveFiles));
% 
% %get the reference electrograms
% thisFilename = 'Wave_refs.csv';
% [refInfo, refVarnames, refData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
% iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
% loadedFiles(iThisFile) = true;
% 
% 
% %get the roving bipolar electrograms
% thisFilename = 'Wave_rov.csv';
% [rovInfo, rovVarnames, rovData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
% iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
% loadedFiles(iThisFile) = true;
% 
% %import the unipolar electrograms
% switch type
%     case 'standard'
%         %import uni distal
%         thisFilename = 'Wave_uni_distal.csv'
%         [uniDistInfo, uniDistVarnames, uniDistData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
%         iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
%         loadedFiles(iThisFile) = true;
% 
%         %import uni proximal
%         thisFilename = 'Wave_uni_proximal.csv';
%         [uniProxInfo, uniProxVarnames, uniProxData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
%         iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
%         loadedFiles(iThisFile) = true;
% 
%     case 'omnipolar'
%         %import uni across
%         thisFilename = 'Wave_uni_across.csv';
%         [uniAcrossInfo, uniAcrossVarnames, uniAcrossData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
%         iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
%         loadedFiles(iThisFile) = true;
% 
%         %import uni along
%         thisFilename = 'Wave_uni_along.csv';
%         [uniAlongInfo, uniAlongVarnames, uniAlongData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
%         iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
%         loadedFiles(iThisFile) = true;
% 
%         %import uni corner
%         thisFilename = 'Wave_uni_corner.csv';
%         [uniCornerInfo, uniCornerVarnames, uniCornerData] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() thisFilename]);
%         iThisFile = find(~cellfun('isempty', regexp(allFiles, ['^' thisFilename], 'once')));
%         loadedFiles(iThisFile) = true;
% end

% %check for any other wave or map files
% if any(~loadedFiles)
%     % we have additional wave files, check whether to just load these or
%     % ask the user what to do
%     if loadallwavefiles
%         extraFilesToLoad = allFiles(~loadedFiles);
%         for iFile = 1:numel(extraFilesToLoad)
%             if strcmpi(extraFilesToLoad{iFile}, 'Wave_refs2.csv')
%                 % This is to make sure that we do not attempt to load a
%                 % Wave_refs2.csv file, which for now seemt to be empty.
%                 % TODO: revisit Wave_refs2.csv files in the future if new
%                 % data is present
%                 continue
%             else
%                 [extraFilesInfo{iFile}, extraFilesVarnames{iFile}, extraFilesData{iFile}] = loadensitex_dxldata([latMapDir{:} filesep() 'Contact_Mapping' filesep() extraFilesToLoad{iFile}]);
%             end
%         end
%     else
%         warning('IMPORTENSITEX_OPENEP: Extra wave files identified, but code to ask the user what to do has not yet been implemented. For now if you want access to these wavefiles, re-run this programme with the option loadallwavefiles set to TRUE');
%     end
% end








%% Calculate annotation times

% first we need to find out which file stored in mappingData is labelled as
% a local activation time map.
isLAT = cellfun(@(s) contains(s.info.mapType, 'LAT'), mappingData);
if sum(isLAT)>1
    error(['IMPORTENSITEX_OPENEP: Too many local activation time mapping files found in folder ' wavesFolder '. Please ensure only one Map_LAT_*.csv file is present.']);
end

% these are all in samples
refTick_adj   = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'Ref Tick'))); % _adj because these already reflect the user adjustments
rovTick_adj   = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'Rov Tick 1')));

%TODO: CHECK THESE TIMES ARE ADJUSTED APPROPRIATELY
leftCurtain     = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'left curtain (ms)'))) / 1000 * rovInfo.sampleFreq;
rightCurtain    = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'right curtain (ms)'))) / 1000 * rovInfo.sampleFreq;

startTime_s       = str2double(rovData(:,strcmpi(rovVarnames, 'startTime (abs)'))); % this comes from the roving wave file
refTime_s         = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'refTime (abs)'))); % this comes from the mapping file
adjTime_ms        = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'adjTime (ms)'))); % from the mapping file
lat_ms            = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'LAT'))); % from the mapping file
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

% Calculate the windows of interest
leftCurtain = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'left curtain (ms)'))) / 1000 * rovInfo.sampleFreq;
rightCurtain = str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'right curtain (ms)'))) / 1000 * rovInfo.sampleFreq;








%% Calculate voltages

isPP = cellfun(@(s) contains(s.info.mapType, 'PP'), mappingData);
if sum(isPP)>1
    error(['IMPORTENSITEX_OPENEP: Too many peak to peak mapping files found in folder ' wavesFolder '. Please ensure only one Map_PP_*.csv file is present.']);
end

bipolarVoltages = str2double(mappingData{isPP}.data(:,strcmpi(mappingData{isPP}.varnames, 'P-P')));
includeFlag = str2double(mappingData{isPP}.data(:,strcmpi(mappingData{isPP}.varnames, 'utilized')));
pointNumberFromFile = mappingData{isPP}.data(:,strcmpi(mappingData{isPP}.varnames, '(Point #)'));








%% Save all data in the OpenEP format

% General data
userdata = openep_createuserdata();
userdata.systemName = 'ensitex';
userdata.notes{1} = [date() ': Created'];

% this is the directory containing the Contact_Mapping folder that the 
% electrograms came from; noting that Contact_Mapping_Model.xml files might
% have been parsed from adjacent directories.
userdata.notes{end+1} = [date() ': userdata.ensitexFolder stores the directory containing the Contact_Mapping folder that was parsed. Additional Contact_Mapping_Model.xml files might have been parsed from adjacent directories.'];
userdata.ensitexFolder = fileparts(T(mapID,:).egmfiles{egmID}); 
userdata.electric.sampleFrequency = rovInfo.sampleFreq;

% Geometry
userdata.surface.triRep = t;
userdata.surface.normals = normals;

surfaceOfOrigin = data_geometry.dxgeo.surface_of_origin;
userdata = setSurfaceProperty(userdata, 'name', 'surfaceOfOrigin', 'map', surfaceOfOrigin, 'definedOn', 'elements');


% Surface maps, removing invalid data beyond interpolation distance
userdata.surface.act_bip = [act bip];
userdata.surface.uni_imp_frc = [uni imp frc];

% Electric data - this should be refactored and moved higher in the code.
% In this section we should only have pre-calculated variables and be
% storing them in userdata, for clarity.
userdata.electric.electrodeNames_bip    = mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'Rov trace'));
userdata.electric.egmX                  = [str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'roving x'))) ...
    str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'roving y'))) ...
    str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'roving z')))];
userdata.electric.egmSurfX              = [str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'surface x'))) ...
    str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'surface y'))) ...
    str2double(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames, 'surface z')))];
userdata.electric.egmRef                = local_concatdata(refData(:,strcmpi(refVarnames,'signals')) ...
    ,refData(:,strcmpi(refVarnames,'Freeze Grp #')) ...
    ,rovData(:,strcmpi(rovVarnames,'Freeze Grp #')) ...
    ,refInfo.filename);
userdata.electric.egm                   = local_concatdata(rovData(:,strcmpi(rovVarnames,'signals')),[],[],rovInfo.filename);

userdata.electric.annotations.referenceAnnot = refAnnot;
userdata.electric.annotations.mapAnnot  = latAnnot;
userdata.electric.annotations.woi       = leftCurtain;
userdata.electric.annotations.woi(:,2)  = rightCurtain;

%userdata.electric.annotations.woi       = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'left curtain (ms)'))) / 1000 * userdata.electric.sampleFrequency;
%userdata.electric.annotations.woi(:,2)  = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'right curtain (ms)'))) / 1000 * userdata.electric.sampleFrequency;

userdata.electric.voltages.bipolar   = bipolarVoltages;
userdata.electric.include            = includeFlag;
userdata.electric.names               = pointNumberFromFile;

% userdata.electric.voltages.bipolar      = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'P-P')));
% userdata.electric.include               = str2double(mappingPoints.data(:,strcmpi(mappingPoints.varnames, 'utilized')));
% userdata.electric.names                 = mappingPoints.data(:,strcmpi(mappingPoints.varnames, '(Point #)'));

% Unipolar electrograms
userdata.electric.electrodeNames_uni    = local_parseuninames(mappingData{isLAT}.data(:,strcmpi(mappingData{isLAT}.varnames,'Electrodes')));
userdata.electric.egmUni = zeros([size(userdata.electric.egm) size(userdata.electric.egmUniX,3)]);

% Unipolar electrogram locations are stored differently for standard and omnipolar configurations
switch egmtype
    case 'bi'
        disp('IMPORTENSITEX_OPENEP: Parsing unipolar co-ordinates for bipolar configuration ...');

        warning('IMPORTENSITEX_OPENEP: When a map is exported we are not given the individual unipole co-ordinates ... assuming uni distal and uni proximal are at the same location')
        userdata.electric.egmUniX = cat(3, userdata.electric.egmX, userdata.electric.egmX);

        disp('IMPORTENSITEX_OPENEP: Parsing unipolar electrograms for bipolar configuration ...');
        userdata.electric.egmUni(:,:,1) = local_concatdata(uniDistData(:,strcmpi(uniDistVarnames,'signals')),[],[],uniDistInfo.filename);
        userdata.electric.egmUni(:,:,2) = local_concatdata(uniProxData(:,strcmpi(uniProxVarnames,'signals')),[],[],uniProxInfo.filename);

    case 'omni'
        disp('IMPORTENSITEX_OPENEP: Parsing unipolar co-ordinates for omnipolar configuration ...');

        % now we need to locate the correct co-ordinates from mappingPoints,
        % without assuming that the order is correct. The order is however
        % usually (1) corner, (2) along, (3) across, so that is what we will
        % check first
        uni1Assigned = false;
        uni2Assigned = false;
        uni3Assigned = false;
        for iPoint = 1:size(userdata.electric.electrodeNames_uni,1)
            if strcmpi(userdata.electric.electrodeNames_uni{iPoint,1}, mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_Corner_Elec')))
                X = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_CornerX'));
                Y = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_CornerY'));
                Z = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_CornerZ'));
                userdata.electric.egmUniX(iPoint,1:3,1) = str2double([X Y Z]);

                % Identify the correct uni data and varnames - corner
                if ~uni1Assigned
                    uni1Data = uniCornerData;
                    uni1Varnames = uniCornerVarnames;
                    uni1Info = uniCornerInfo;
                    uni1Assigned = true;
                end
            else
                error('IMPORTENSITEX_OPENEP: Electrode naming mismatch')
            end
            if strcmpi(userdata.electric.electrodeNames_uni{iPoint,2}, mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_Along_Elec')))
                X = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_AlongX'));
                Y = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_AlongY'));
                Z = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_AlongZ'));
                userdata.electric.egmUniX(iPoint,1:3,2) = str2double([X Y Z]);

                % Identify the correct uni data and varnames - along
                if ~uni2Assigned
                    uni2Data = uniAlongData;
                    uni2Varnames = uniAlongVarnames;
                    uni2Info = uniAlongInfo;
                    uni2Assigned = true;
                end
            else
                error('IMPORTENSITEX_OPENEP: Electrode naming mismatch')
            end
            if strcmpi(userdata.electric.electrodeNames_uni{iPoint,3}, mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_Across_Elec')))
                X = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_AcrossX'));
                Y = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_AcrossY'));
                Z = mappingData{isLAT}.data(iPoint,strcmpi(mappingData{isLAT}.varnames,'Uni_AcrossZ'));
                userdata.electric.egmUniX(iPoint,1:3,3) = str2double([X Y Z]);

                % Identify the corret uni data and varnames - across
                if ~uni3Assigned
                    uni3Data = uniAcrossData;
                    uni3Varnames = uniAcrossVarnames;
                    uni3Info = uniAcrossInfo;
                    uni3Assigned = true;
                end
            else
                error('IMPORTENSITEX_OPENEP: Electrode naming mismatch')
            end
        end

        userdata.electric.egmUni(:,:,1) = local_concatdata(uni1Data(:,strcmpi(uni1Varnames,'signals')),[],[],uni1Info.filename);
        userdata.electric.egmUni(:,:,2) = local_concatdata(uni2Data(:,strcmpi(uni2Varnames,'signals')),[],[],uni2Info.filename);
        userdata.electric.egmUni(:,:,3) = local_concatdata(uni3Data(:,strcmpi(uni3Varnames,'signals')),[],[],uni3Info.filename);

        %userdata.electric.egmUniSurfX           = userdata.electric.egmUniX; % note that we do not have co-ordinates for the second unipole
        disp('IMPORTENSITEX_OPENEP: Finished parsing unipolar co-ordinates ...');
end

% % Store any additional signals in the ecg array
% disp('dealing with ECG electrograms')
% userdata.electric.ecgNames = {};
% if ~isempty(extraFilesInfo)
%     % for speed first work out the dimensions and pre-populate
%     fWait = waitbar(0, 'Storing additional ECG names');
%     for iEF = 1:numel(extraFilesInfo)
%         if any(strcmpi(extraFilesVarnames{iEF}, 'signals'))
%             % Then, this extra file contains signal data - store these in
%             % the ECG array.
% 
%             % First check for unique electrode names in this file
%             electrodes = unique(extraFilesData{iEF}(:,1));
% 
%             % Remove any non-ASCII characters, leading or trailing spaces
%             % and duplicate rows (cElectrodes for 'clean electrodes')
%             cElectrodes = unique(cellfun(@(s) strtrim(regexprep(s, '[^\x00-\x7F]', '')), electrodes, 'UniformOutput', false));
% 
%             % Add electrode names to the ecgNames cell array
%             userdata.electric.ecgNames = union(userdata.electric.ecgNames, cElectrodes);
%         else
%             % Then, this extra file does not contain signal data. If the
%             % file has not already been imported (we do not yet have a
%             % check for this) then it is likely to be an additional mapping
%             % file. Do nothing for the time being.
%         end
%         waitbar(iEF/numel(extraFilesInfo),fWait);
%     end
%     close(fWait);
% 
%     % prepopulate for speed
%     userdata.electric.ecg = zeros([size(userdata.electric.egm) size(userdata.electric.ecgNames,1)]);
% 
%     fWait = waitbar(0, 'Storing additional ECG data');
%     for iEF = 1:numel(extraFilesInfo)
%         % check if this is a signals file
%         if any(strcmpi(extraFilesVarnames{iEF}, 'signals'))
%             % Next we iterate through every signal and work out where to
%             % put it in the ECG array.
% 
%             for jSg = 1:size(extraFilesData{iEF},1)
%                 thisSig = extraFilesData{iEF}(jSg,strcmpi(extraFilesVarnames{iEF},'signals'));
%                 thisName = extraFilesData{iEF}(jSg,strcmpi(extraFilesVarnames{iEF},'Trace'));
% 
%                 % clean the name
%                 thisName = strtrim(regexprep(thisName, '[^\x00-\x7F]', ''));
% 
%                 userdata.electric.ecg(jSg,:,strcmpi(userdata.electric.ecgNames, thisName)) = thisSig{:};
%             end
%         else
%             % Then, this extra file does not contain signal data. If the
%             % file has not already been imported (we do not yet have a
%             % check for this) then it is likely to be an additional mapping
%             % file. Do nothing for the time being.
%         end
%         waitbar(iEF/numel(extraFilesInfo),fWait);
%     end
%     close(fWait);
% end

% set up the surface normals
tr = getMesh(userdata, 'triangulation');
[closestVertices,~] = findclosestvertex(tr, userdata.electric.egmX, true);
userdata.electric.barDirection = userdata.surface.normals(closestVertices,:);

% we don't have impedance values, so create NaN values
userdata.electric.impedances.time = cell(7110,1);
userdata.electric.impedances.value = cell(7110,1);
[userdata.electric.impedances.value{:}] = deal(NaN);
[userdata.electric.impedances.time{:}] = deal(NaN);

% we don't have the unipolar peak to peak voltages so we have to do something
userdata.electric.voltages.unipolar = NaN(size(userdata.electric.voltages.bipolar));
userdata.electric.voltages.unipolar = calculatePeak2PeakVoltage( userdata.electric.egmUni, userdata.electric.annotations.referenceAnnot, userdata.electric.annotations.woi );

% Temp - remove signalMaps which, if empty, prevents the file being loaded in EP Workbench
userdata.surface = rmfield(userdata.surface, 'signalMaps');
userdata.electric.tags = cell(length(userdata.electric.names),1);










%% Encourage user to save the data
matFileFullPath = [];
if ~isempty(saveFileName)
    save(saveFileName, 'userdata');
    matFileFullPath = saveFileName;
else
    defaultName = [mappingData{isLAT}.info.study '_' mapToRead];
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
        tf = cellfun(@(p) contains(p, stub, 'IgnoreCase', true), allSubFolders);
        thisFolder = allSubFolders(tf);
        % Check if more than one folder meets the critiera, and ask the user to choose
        if numel(thisFolder)>1
            warning(['IMPORTENSITEX_OPENEP: More than one candidate folder selected for the export of ***' stub '*** data. Please choose one folder ...'])
            [indx, tf] = listdlg('ListString', thisFolder ...
                ,'ListSize', [480 300] ...
                , 'name', ['Which is the correct ***' stub '*** folder?'] ...
                , 'selectionmode', 'single' ...
                );
            if ~tf
                error('IMPORTENSITEX_OPENEP: Operation cancelled')
            else
                thisFolder = thisFolder{indx};
            end
        end
        pathName = fullfile(studyDir, thisFolder);
        pathName = pathName{:};
    end

    function hd = local_homedirec()
        %HOMEDIREC returns the user's home directory.

        if ispc
            hd = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
        else
            hd = getenv('HOME');
        end
    end

    function matrixData = local_concatdata(cellData, freezeGroupIn, freezeGroupOut, dataName)
        % This function concatenates cell data into a matrix, optionally
        % based on the ordering specified by freezeGroupIn and freezeGroupOut
        f = waitbar(0, ['Reorganising data for:' dataName]);
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

        % pre-allocate for speed
        matrixData = zeros(nCell,size(cellDataNew{1},2)); % we assume that all cells have the same length
        % matrixData = cellDataNew{1};
        for iCell = 1:nCell
            matrixData(iCell,:) = cellDataNew{iCell};
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
                if ~isempty(splt{j})
                    uniNames{iPair,j} = splt{j};
                end
            end
            % uniNames{iPair,1} = splt{1};
            % uniNames{iPair,2} = splt{2};
        end
    end

    function xmlFiles = local_findAllXmlFiles(parentDirectory)
        % local_findAllXmlFiles  Recursively finds all .xml files under parentDirectory.
        
        % Use dir with recursive wildcard
        fileList = dir(fullfile(parentDirectory, '**', '*.xml'));

        % Filter hidden files
        fileList = fileList(~startsWith({fileList.name}, '.'));

        % Extract full paths into a cell array
        xmlFiles = fullfile({fileList.folder}, {fileList.name});
    end

    function csvFiles = local_findAllCsvFiles(parentDirectory)
        % local_findAllCsvFiles  Recursively finds all .csv files under parentDirectory

        % Use dir with recursive wildcard
        fileList = dir(fullfile(parentDirectory, '**', '*.csv'));

        % Filter hidden files
        fileList = fileList(~startsWith({fileList.name}, '.'));

        % Extract full paths into a cell array
        csvFiles = fullfile({fileList.folder}, {fileList.name});
    end

    function info = local_loadCsvFileHeader(csvFile)

        info = [];
        fileID = fopen(csvFile, 'r');
        if fileID == (-1)
            error('LOADENSITEX_DXLDATA: Could not open file.')
        end
        cleanupFile = onCleanup(@()fclose(fileID));

        maxBytes = 100000; % enough data to cover the header
        fseek(fileID, 0, 'bof'); % move to the beginning of the file
        if maxBytes > filebytes2end(fileID)
            maxBytes = filebytes2end(fileID);
        end
        [fData, fDataSize] = fread(fileID, maxBytes, '*char');
        fData = fData(1:fDataSize)';

        % do the prechecks and return if bad
        if ~loadensitex_prechecks(fData, 'DxL')
            warning('IMPORTENSITEX_OPENEP/LOCAL_LOADCSVFILEHEADER: A non-DxL CSV file was identified, such as a lesions or automark file. We will ignore this.')
            return
        end

        % The 'header' finishes at the end of the last line starting with "****,"
        [~, ind2] = regexp(fData, '****','start','end');
        if isempty(ind2)
            error('End of header not found. Double check that maxBytes is large enough to cover header.')
        end
        indEndofHeader = ind2(end);
        header = fData(1:indEndofHeader);

        % Parse the header
        info = parse_header(header, 'dxl');
        info.filename = csvFile;

    end

    function tf = local_compareXmlFiles(S1, S2)
        % compareMeshStructs  Compare two mesh structures containing dxgeo subfields.
        %
        % Returns true only if S1.dxgeo and S2.dxgeo both contain the fields
        % 'vertices', 'triangles', and 'normals', and all three arrays are exactly equal.

        % Required subfields within dxgeo
        requiredFields = {'vertices', 'triangles', 'normals'};

        % Check dxgeo exists in both structures
        if ~isfield(S1, 'dxgeo') || ~isfield(S2, 'dxgeo')
            tf = false;
            return;
        end

        % Check required subfields exist
        for k = 1:numel(requiredFields)
            f = requiredFields{k};
            if ~isfield(S1.dxgeo, f) || ~isfield(S2.dxgeo, f)
                tf = false;
                return;
            end
        end

        % Compare arrays for exact equality
        tf = isequal(S1.dxgeo.vertices,  S2.dxgeo.vertices)  && ...
            isequal(S1.dxgeo.triangles, S2.dxgeo.triangles) && ...
            isequal(S1.dxgeo.normals,   S2.dxgeo.normals);
    end

    function [tf, identifiedType] = local_areMapSuffixesUniform(allMapTypes)
        % Return true if all map suffixes (text after last '_')
        % are identical in all entries of the input cell array of strings.

        % Extract suffix from each map type
        suffixes = cellfun( ...
            @(s) s( find(s=='_',1,'last')+1 : end ), ...
            allMapTypes, 'UniformOutput', false);

        % True only if all suffixes are the same
        tf = numel(unique(suffixes)) == 1;

        % The type, if true
        if tf
            identifiedType = suffixes{1};
        else
            identifiedType = false;
        end
    end

    function [tf, numMappingPts] = local_areMapPointNumbersUniform(allMaps)
        % Return true if all maps in the input cell array of map headers 
        % have the same number of mapping points. Return these in numMappingPoints

        % Extract the number of mapping points
        numPoints = cellfun(@(s) s.numPoints, allMaps);

        % True only if all numbers are the same
        tf = numel(unique(numPoints)) == 1;
        
        if tf
            numMappingPts = numPoints(1);
        else
            numMappingPts = false;
        end
    end

    function out = local_lastTwoParts(paths)
        % lastTwoParts returns the last two parts of each path in a cell array
        % paths: 1×N or N×1 cell array of full paths
        % out:   1×N cell array of "secondLastPart/lastPart"

        if ischar(paths) || isstring(paths)
            paths = {char(paths)};
        elseif ~iscell(paths)
            error('Input must be char, string, or cell array of char/string.');
        end

        N = numel(paths);
        out = cell(1, N);

        for k = 1:N
            p = paths{k};

            % Get the last part
            [~, lastPart, ext] = fileparts(p);
            lastPartFull = [lastPart ext];  % include extension if any

            % Get the second-to-last folder
            [parentFolder, secondLastPart] = fileparts(fileparts(p));

            % Combine
            out{k} = fullfile(secondLastPart, lastPartFull);
        end
    end



 

end