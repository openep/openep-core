function [userdata, matFileFullPath] = importcarto_mem(varargin)
% IMPORTCARTO provides a data structure from multiple carto files (from zip).
% Usage:
%   userdata = importcarto_mem(userinput)
%   userdata = importcarto_mem()
%   [userdata, matFileFullPath] = ...
% Where:
%   dirName is the directory with all of the files corresponding to a map
%   userdata is a single data structure
%   matFileFullPath is the path to the .mat file, if opened or saved
%
% IMPORTCARTO can load data in the following ways:
%   *) USERINPUT is a .mat file containing userdata
%   *) USERINPUT is a .xml file - this must be the xml file in a folder
%                                   containing all the other Carto3 files.
%
% IMPORTCARTO_MEM accepts the following parameter-value pairs
%   'maptoread'     {''}|string|double
%       Specifies which map to read. Can be a string referring
%       to the map name or a double referring to the number of points in the
%       map. If there are multiple maps with the same number of points an error
%       will be thrown.
%   'refchannel'    {''}|string
%       The name of the channel to pick as the refence channel. Typically
%       this is the pacing channel for the map. Specify a string such as
%       'CS9-CS10'.
%   'ecgchannel'    {''}|string
%       The name of the channel to pick as the ECG channel. Typically
%       this is an informative ECG such as V1. Specify a string such as
%       'V1'.
%   'savefilename'       {''}|string
%       The full path to the location in which to save the output.
%   'verbose'       {true}|false
%       Not yet implemented
% Example of command line entry ...
%       userdata = importcarto_mem(<path to XML file>, ...
%                                        'maptoread', 1693,  ...
%                                        'refchannel', 'CS9-CS10', ...
%                                        'ecgchannel', 'V1')
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

% Author: Nick Linton (2011)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications - Steven Williams (2013) - force data
% Modifications - Steven Williams (2018) - unipole data
% Modifications - Steven Williams (2020) - command line input
% Modifications - Nick Linton (2023)
%   - cleaned, removed need to check the position of electrodes
%   - changed unipolar data to be just the one pole referenced in the
%     ECG_Export file
%
% ---------------------------------------------------------------
%

% testing
% hSurf = drawMap(userdata, 'type', 'bip', 'coloraxis', [0 2])
% set(hSurf, 'facealpha', 0.5)
% x = userdata.electric.egmX(:,1);
% y = userdata.electric.egmX(:,2);
% z = userdata.electric.egmX(:,3);
% x1 = userdata.electric.egmUniX(:,1,2);
% y1 = userdata.electric.egmUniX(:,2,2);
% z1 = userdata.electric.egmUniX(:,3,2);
% hold on
% hX1 = plot3(x,y,z,'.','markersize', 10, 'color', 'k')
% hX2 = plot3(x1,y1,z1,'.','markersize', 10, 'color', 'g')
% ---------------------------------------------------------------

% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

persistent saveDir homeDir

matFileFullPath = [];

if isempty(saveDir) || ~ischar(saveDir) || ~isfolder(saveDir)
    if exist('stevematlabroot.m','file')==2
        saveDir = stevematlabroot();
    elseif exist('nickmatlabroot.m','file')==2
        saveDir = nickmatlabroot();
    else
        saveDir = matlabroot();
    end
end

if isempty(homeDir) || ~isfolder(homeDir);     homeDir = saveDir; end

userdata = [];
hWait = [];

if nargin >= 1
    userinput = varargin{1};
else
    dialog_title = 'Select the Carto Study xml in the unzipped folder (eg "Study 1 11_20_2012 21-02-32.xml"), or mat file.';
    filterSpec = {'*.zip;*.xml;*.mat', 'Appropriate files (*.zip;*.xml;*.mat)' ; '*.zip','Zip (*.zip)' ; '*.xml','XML (*.xml)' ; '*.mat','Matlab (*.mat)' ; '*.*','All files (*.*)'};
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
    error('IMPORTCARTO_MEM: ZIP files are no longer supported - unzip externally to matlab.');
    return %#ok<UNRCH>
elseif strcmpi(userinput((end-3):end),'.mat')
    s = load(userinput);
    userdata = s.userdata;
    matFileFullPath = userinput;
    return
else
    studyDir = fileparts(userinput);
    homeDir = studyDir;
end

% parse command line input
nStandardArgs = 1; % UPDATE VALUE
mapToRead_cli = '';
channelRef_cli = '';
channelECG_cli = '';
saveFileName_cli = '';
verbose = true;
if nargin > nStandardArgs
    for i = nStandardArgs+1:2:nargin
        switch varargin{i}
            case 'maptoread'
                mapToRead_cli = varargin{i+1};
            case 'refchannel'
                channelRef_cli = varargin{i+1};
            case 'ecgchannel'
                channelECG_cli = varargin{i+1};
                if ischar(channelECG_cli); channelECG_cli = {channelECG_cli}; end
            case 'savefilename'
                saveFileName_cli = varargin{i+1};
            case 'verbose'
                verbose = varargin{i+1};
            otherwise
                error('IMPORTCARTO_MEM: Unrecognized input.')
        end
    end
end

% try
studyDirInfo = dir(studyDir);
allfilenames = cell(length(studyDirInfo),1);
for i = 1:length(studyDirInfo)
    allfilenames{i} = studyDirInfo(i).name;
end

disp(['Accessing: ' userinput]);

hWait = waitbar(0, 'Getting study information');

Pref.Str2Num = 'never';
[tree, ~, ~] = xml_read(userinput, Pref);

delete(hWait)
hWait = [];

study = tree.ATTRIBUTE.name; %#ok<NASGU>
nMaps = numel(tree.Maps.Map);

names = cell(nMaps, 1);
for iMap = 1:nMaps
    names{iMap} = tree.Maps.Map(iMap).ATTRIBUTE.Name;
    numPtsPerMap(iMap) = str2double(tree.Maps.Map(iMap).CartoPoints.ATTRIBUTE.Count);
end

if isempty(mapToRead_cli)
    [selection,ok] = listdlg(     'ListString', names ...
        , 'SelectionMode', 'single' ...
        , 'PromptString', 'Which map do you want to access?' ...
        , 'ListSize', [300 300] ...
        );
    if ~ok
        return
    end
else
    if isnumeric(mapToRead_cli)
        selection = find(numPtsPerMap==mapToRead_cli);
        if numel(selection)>1
            error(['IMPORTCARTO_MEM: Multiple maps with ' ...
                num2str(mapToRead_cli) ...
                ' points identified. Use an alternative method to identify map.']);
        elseif isempty(selection)
            error(['IMPORTCARTO_MEM: No map with ' ...
                num2str(mapToRead_cli) ...
                ' points identified. Check the number of points specified is correct.']);
        end
    elseif ischar(mapToRead_cli)
        selection = find(strstartcmpi(mapToRead_cli, names));
    end
end

% Get the tags from the ID
nTags = str2double(tree.Maps.TagsTable.ATTRIBUTE.Count);
tagNames = cell(nTags,1);
tagID = zeros(nTags,1);
for iTag = 1:nTags
    tagNames{iTag} = tree.Maps.TagsTable.Tag(iTag).ATTRIBUTE.Full_Name;
    tagID(iTag) = str2double(tree.Maps.TagsTable.Tag(iTag).ATTRIBUTE.ID);
end

map = [];
for iMap = selection
    cartoMap = tree.Maps.Map(iMap);
    map.studyName = tree.ATTRIBUTE.name;
    map.name = cartoMap.ATTRIBUTE.Name;
    map.meshFile = cartoMap.ATTRIBUTE.FileNames;


    %%% Work out the root filename
    filenameroot = map.meshFile;
    k = strfind(filenameroot, '.mesh');
    filenameroot(k:end) = [];


    isFirstPointRead = false;
    nPoints = str2double(cartoMap.CartoPoints.ATTRIBUTE.Count);
    if nPoints>0
        filename = [filenameroot '_P' num2str(cartoMap.CartoPoints.Point(1).ATTRIBUTE.Id) '_ECG_Export.txt'];
        filename = mycheckfilename(filename, allfilenames, ['P' num2str(cartoMap.CartoPoints.Point(1).ATTRIBUTE.Id) '_ECG_Export']);
        if ~isempty(filename)
            ecgFileHeader = read_ecgfile_v4(fullfile(studyDir, filename));
            names = ecgFileHeader.channelNames;
            namesFull = ecgFileHeader.channelNamesFull;

            if isempty(channelRef_cli)
                [kRef,ok] = listdlg( 'ListString', names , 'SelectionMode','single' , 'PromptString','Which signal is Ref?' , 'ListSize',[300 300] ); if ~ok; return; end
                channelRef_cli = names{kRef};
            else
                kRef = find(strstartcmpi(channelRef_cli, names));
                if isempty(kRef) || numel(kRef)>1
                    error(['IMPORTCARTO_MEM: Unable to uniquely identify the specified reference channel: ' channelRef_cli]);
                end
            end
            if isempty(channelECG_cli)
                [kEcg,ok] = listdlg( 'ListString', names , 'SelectionMode','multiple' , 'PromptString','Which other signals should be downloaded with each point (typically one or more ECG signals)?' , 'ListSize',[300 300] ); if ~ok; return; end
                channelECG_cli = names(kEcg);
            else
                kEcg = zeros(1,numel(channelECG_cli));
                for i = 1:numel(channelECG_cli)
                    kEcg(i) = find(strstartcmpi(channelECG_cli{i}, names));
                    if isempty(kEcg(i))
                        error(['IMPORTCARTO_MEM: Unable to uniquely identify the specified ECG channel: ' channelECG_cli]);
                    end
                end
            end
            isFirstPointRead = true;
        end
    end




    nAnat = 0;
    if isfield(cartoMap, 'Anatomical_Tags')
        if ~isempty(cartoMap.Anatomical_Tags)
            for iAnat = 1:length(cartoMap.Anatomical_Tags.Anatomical_Tag)
                if isempty(cartoMap.Anatomical_Tags.Anatomical_Tag(iAnat).Perimiter.LINE_LOOP.CONTENT)
                    continue;
                else
                    temp = str2num(cartoMap.Anatomical_Tags.Anatomical_Tag(iAnat).Perimiter.LINE_LOOP.CONTENT);
                end
                if ~isempty(temp)
                    nAnat = nAnat + 1;
                    map.anatPerimeter{nAnat} =  reshape(temp,3,numel(temp)/3)';
                end
            end
        end
    end


    % Now get the new points in this map
    nPoints = str2double(cartoMap.CartoPoints.ATTRIBUTE.Count);
    egm = [];
    egmUni1 = [];
    egmUni2 = [];
    ref = [];
    ecg = [];
    map.xyz = NaN(nPoints, 3);
    map.xyz2 = NaN(nPoints, 3);
    %map.xyzSurf = NaN(nPoints, 3);  %not reliable
    %map.projDist = NaN(nPoints, 3); %not reliable
    map.pointNames = cell(nPoints, 1);
    map.pointID = zeros(nPoints, 1);
    map.tag = cell(nPoints, 1);
    map.electrodesUsed = cell(nPoints, 1);


    if nPoints>0
        for iPoint = 1:nPoints
            map.xyz(iPoint,:) = str2num(cartoMap.CartoPoints.Point(iPoint).ATTRIBUTE.Position3D);
            %map.xyzSurf(iPoint,:) = cartoMap.CartoPoints.Point(iPoint).VirtualPoint.ATTRIBUTE.Position3D;   %not reliable
            %map.projDist(iPoint,:) = cartoMap.CartoPoints.Point(iPoint).VirtualPoint.ATTRIBUTE.ProjectionDistance;  %not reliable
            map.pointID(iPoint) = str2double(cartoMap.CartoPoints.Point(iPoint).ATTRIBUTE.Id);
            map.pointNames{iPoint} = [  'P' num2str(map.pointID(iPoint))  ];
            temp = [];
            if isfield(cartoMap.CartoPoints.Point(iPoint), 'Tags')
                temp = cartoMap.CartoPoints.Point(iPoint).Tags;
            end
            if ~isempty(temp)
                content = str2num(temp.CONTENT);
                iTag = find(content==tagID,1,'first');
                if numel(content)>1
                    beep
                    warning('Nick to reprogram for multiple tags')
                end
                map.tag{iPoint} = tagNames{iTag};
            end
        end
    end

    %%% Read in the ContactForceInRF_ and RF_ files (if they exist)
    iRfFiles = find(strstartcmpi(['RF_' map.name], allfilenames));
    iContactForceInRfFiles = find(strstartcmpi(['ContactForceInRF_' map.name], allfilenames));
    C=0;R=1;
    if ~isempty(iRfFiles)
        for i = 1:numel(iRfFiles)
            if i == 1
                rfData = dlmread(fullfile(studyDir, allfilenames{iRfFiles(i)}),'\t',R,C); %first RF_ file
            else
                tempdata = dlmread(fullfile(studyDir, allfilenames{iRfFiles(i)}),'\t',R,C);
                rfData(end+1:end+size(tempdata,1),:) = tempdata; %subsequent RF_ files
            end
        end
        rfData(:,end+1) = NaN; %preallocate the column that will be used for indexing position
        numColRfData = size(rfData,2);
    end
    if ~isempty(iContactForceInRfFiles)
        for i = 1:numel(iContactForceInRfFiles)
            if i == 1
                contactForceInRfData = dlmread(fullfile(studyDir, allfilenames{iContactForceInRfFiles(i)}),'\t',R,C); %first ContactForceInRf_ file
            else
                tempdata = dlmread(fullfile(studyDir, allfilenames{iContactForceInRfFiles(i)}),'\t',R,C);
                contactForceInRfData(end+1:end+size(tempdata,1),:) = tempdata; %subsequent ContactForceInRf_ files
            end
        end
        contactForceInRfData(:,end+1) = NaN; %preallocate the column that will be used for indexing position
        numColCfInRfData = size(contactForceInRfData,2);
    end

    %%% Now get the point WOI, Reference time and Annotation time
    hWait = waitbar(0, ['Getting annotation data for ' num2str(nPoints) ' points']);
    pointExport_WOI = NaN(nPoints,2);
    pointExport_ReferenceAnnotation = NaN(nPoints,1);
    pointExport_MapAnnotation = NaN(nPoints,1);
    pointExport_Unipolar = NaN(nPoints,1);
    pointExport_Bipolar = NaN(nPoints,1);
    electrodeNames_bip = cell(nPoints, 1);
    electrodeNames_uni = cell(nPoints, 2);
    unipolarEgms = NaN(size(egmUni1,1), size(egmUni1,2), 2);
    unipolarEgmsX = NaN(size(map.xyz,1), size(map.xyz,2), 2);
    pointExport_ImpedanceTime = cell(nPoints,1);
    pointExport_ImpedanceValue = cell(nPoints,1);

    if nPoints>0
        for iPoint = 1:nPoints
            waitbar(iPoint/nPoints, hWait);
            filename_pointExport = [filenameroot '_' map.pointNames{iPoint} '_Point_Export.xml'];
            if ~isfile(fullfile(studyDir, filename_pointExport))
                disp(['File not found: ' , filename_pointExport])
            else
                [pointExportTree, ~, ~] = xml_read(fullfile(studyDir, filename_pointExport), Pref);
                pointExport_WOI(iPoint,:) = [str2double(pointExportTree.WOI.ATTRIBUTE.From) str2double(pointExportTree.WOI.ATTRIBUTE.To)];
                pointExport_ReferenceAnnotation(iPoint) = str2double(pointExportTree.Annotations.ATTRIBUTE.Reference_Annotation);
                pointExport_MapAnnotation(iPoint) = str2double(pointExportTree.Annotations.ATTRIBUTE.Map_Annotation);
                pointExport_Unipolar(iPoint) = str2double(pointExportTree.Voltages.ATTRIBUTE.Unipolar);
                pointExport_Bipolar(iPoint) = str2double(pointExportTree.Voltages.ATTRIBUTE.Bipolar);

                %TODO: Put an if statement here for if number=0
                try
                    pointExport_ImpedanceTime{iPoint} = arrayfun(@(x) str2double(x.ATTRIBUTE.Time), pointExportTree.Impedances.Impedance);
                    pointExport_ImpedanceValue{iPoint} = arrayfun(@(x) str2double(x.ATTRIBUTE.Value), pointExportTree.Impedances.Impedance);
                catch
                    pointExport_ImpedanceTime{iPoint} = NaN;
                    pointExport_ImpedanceValue{iPoint} = NaN;
                end
            end
        end
        delete(hWait)
        hWait = [];

        %%% Now get the details for the xml files of each point.
        [allPointExport, ~, ~] = xml_read(fullfile(homeDir, [map.name '_Points_Export.xml']));

        if nPoints ~= numel(allPointExport.Point)
            error('IMPORTCARTO_MEM: There is a discrepancy in the number of points between files.')
        end

        if isFirstPointRead
            nEgmSamples = ecgFileHeader.nSamples;
            egm = zeros(nPoints, nEgmSamples);
            egmUni1 = zeros(nPoints, nEgmSamples);
            egmUni2 = zeros(nPoints, nEgmSamples);
            ref = zeros(nPoints, nEgmSamples);
            ecg = zeros(nPoints, nEgmSamples, numel(kEcg));

            nameRef = names{kRef};
            nameEcg = names(kEcg);
            nameRefFull = namesFull{kRef};
            nameEcgFull = namesFull(kEcg);

            %%% Now get the electrograms
            hWait = waitbar(0, ['Getting electrical data for ' num2str(nPoints) ' points']);

            for iPoint = 1:nPoints
                waitbar(iPoint/nPoints , hWait )
                filename = [filenameroot '_' map.pointNames{iPoint} '_ECG_Export.txt'];
                filename = mycheckfilename(filename, allfilenames, [map.pointNames{iPoint} '_ECG_Export']);

                if ~isempty(filename)
                    [headerInfo, voltages] = read_ecgfile_v4(fullfile(studyDir, filename));
                    voltages = voltages * headerInfo.gain;
                    names = headerInfo.channelNames;
                    electrodeNames_bip{iPoint} = headerInfo.bipMapChannel;
                    electrodeNames_uni{iPoint,1} = headerInfo.uniMapChannel;
                    electrodeNames_uni{iPoint,2} = headerInfo.uniMapChannel2;
                    
                    pointFileName = allPointExport.Point(iPoint).ATTRIBUTE.File_Name;
                    pointFileName = fullfile(homeDir, pointFileName);
                    

                    
                    if any(kRef>numel(names)) || any(kEcg>numel(names)) || any(~strcmpi(names(kRef),nameRef)) || any(~strcmpi(names(kEcg),nameEcg))
                        beep()
                        warning(['IMPORTCARTO_MEM: The columns containing data in the .txt files change names. In ' filename])
                        kRef = find( strcmpi(nameRef, names) );
                        if isempty(kRef)
                            warning(['IMPORTCARTO_MEM: The requested reference channel, ' nameRef ' was not found in file: ' filename '. NaN values will be assigned as the reference for this point' ]);
                            kRef = NaN;
                        end
                        for i = 1:numel(kEcg)
                            new_index = find(strstartcmpi(channelECG_cli{i}, names));
                            if isempty(new_index)
                                error('OPENEP:missingChannel', ['IMPORTCARTO_MEM: Unable to uniquely identify the specified ECG channel: ' channelECG_cli{i} 'in file: ' filename]);
                            end
                            kEcg(i) = new_index;
                        end
                        if numel(kMap_bip)~=1 || numel(kRef)~=1 || numel(kEcg)~=numel(nameEcg)
                            error(['IMPORTCARTO_MEM: The name change could not be resolved. Problem with filename: ' filename])
                        end
                    end
                    if ~isempty(electrodeNames_bip{iPoint})
                        kMap_bip = find( strcmpi(electrodeNames_bip{iPoint}, names) );
                        kMap_uni = [find( strcmpi(electrodeNames_uni{iPoint,1}, names) ) find( strcmpi(electrodeNames_uni{iPoint,2}, names) )];
                        egm(iPoint,:) = voltages(:,kMap_bip);
                        egmUni1(iPoint,:) = voltages(:,kMap_uni(1));
                        egmUni2(iPoint,:) = voltages(:,kMap_uni(2));
                        if isnan(kRef)
                            ref(iPoint,:) = NaN;
                        else
                            ref(iPoint,:) = voltages(:,kRef);
                        end
                        ecg(iPoint,:,:) = voltages(:,kEcg);
                        % read in the required positions (for this version
                        % of code, only egm2
                        map.xyz2(iPoint,:) = read_electrodePositionsOnAnnotation(electrodeNames_uni{iPoint,2}, pointFileName);
                    else
                        warning('IMPORTCARTO_MEM: No electrode found ... check "OnAnnotation" file ...')
                        disp(filename)
                        disp('')
                    end
                end
            end
            delete(hWait)
            hWait = [];
        end

        %%% Based on the above data concatenate matrices for the
        %%% unipolar electrograms and their co-ordinates
        unipolarEgms = NaN(size(egmUni1,1), size(egmUni1,2), 2); 
        unipolarEgms(:,:,1) = egmUni1;
        unipolarEgms(:,:,2) = egmUni2;

        unipolarEgmsX = NaN(size(map.xyz,1), size(map.xyz,2), 2); 
        unipolarEgmsX(:,:,1) = map.xyz;
        unipolarEgmsX(:,:,2) = map.xyz2;

        %%% Now get the force data for the first point
        filename_force = [filenameroot '_' map.pointNames{1} '_ContactForce.txt'];
        filename_force = mycheckfilename(filename_force, allfilenames, [map.pointNames{1} '_ContactForce.txt']);
        if ~isempty(filename_force)
            [~, ~, ~, t_T, ~, ~, ~, ~] = read_forcefile_v2(fullfile(studyDir, filename_force));

            force = nan(nPoints,1);
            axialAngle = nan(nPoints,1);
            lateralAngle = nan(nPoints,1);
            t_force = nan(nPoints,max(size(t_T)),2);
            t_axialAngle = nan(nPoints,max(size(t_T)),2);
            t_lateralAngle = nan(nPoints,max(size(t_T)),2);

            %%% Now we get the forces
            hWait = waitbar(0, ['Getting force data for ' num2str(nPoints) ' points']);
            hasWarned = false;
            for iPoint = 1:nPoints
                waitbar(iPoint/nPoints, hWait);
                filename_force = [filenameroot '_' map.pointNames{iPoint} '_ContactForce.txt'];
                filename_force = mycheckfilename(filename_force, allfilenames, [map.pointNames{iPoint} '_ContactForce.txt']);
                if ~isempty(filename_force)
                    [f, aA, lA, t_T, t_F, t_aA, t_lA, systemTime] = read_forcefile_v2(fullfile(studyDir, filename_force));
                    force(iPoint,1) = str2double(f);
                    axialAngle(iPoint,1) = str2double(aA);
                    lateralAngle(iPoint,1) = str2double(lA);
                    t_force(iPoint,1:max(size(t_T)),1) = t_T;
                    t_force(iPoint,1:max(size(t_T)),2) = t_F;
                    t_axialAngle(iPoint,1:max(size(t_T)),1) = t_T;
                    t_axialAngle(iPoint,1:max(size(t_T)),2) = t_aA;
                    t_lateralAngle(iPoint,1:max(size(t_T)),1) = t_T;
                    t_lateralAngle(iPoint,1:max(size(t_T)),2) = t_lA;

                    % add indices to the RF datasets referring back to
                    % egmX and egmSurfX
                    if strstartcmpi('abl', map.tag{iPoint}) && ~isempty(iContactForceInRfFiles)
                        %then we are dealing with an ablation point
                        [~, iMin] = min(abs(t_T));
                        acqTime = systemTime(iMin);

                        if ~isempty(iContactForceInRfFiles) % first for ContactForceInRF data
                            pointSysTime = contactForceInRfData(:,1);
                            [~,iMin] = min(abs(pointSysTime - acqTime));
                            contactForceInRfData(iMin,numColCfInRfData) = iPoint;
                        end

                        if ~isempty(iRfFiles) % then for RF data
                            pointSysTime = rfData(:,1);
                            [~,iMin] = min(pointSysTime - acqTime);
                            rfData(iMin, numColRfData) = iPoint;
                        end
                    end
                elseif ~hasWarned
                    %warning(['IMPORTCARTO_MEM: Force file not found for point ' num2str(iPoint)]);
                    warning(['IMPORTCARTO_MEM: Force file not found for some points' char(10) 'Check userdata.electric.force for missing data.']);
                    hasWarned = true;
                end

            end
            delete(hWait)
            hWait = [];
        end
    else
        nameRef = [];
        nameEcg = [];
        delete(hWait)
        hWait = [];
    end


    % now load into userdata structure
    userdata = [];

    userdata.cartoFolder = studyDir;


    userdata.electric.tags = map.tag;
    userdata.electric.names = map.pointNames;
    userdata.electric.electrodeNames_bip = electrodeNames_bip;
    userdata.electric.egmX = map.xyz;
    userdata.electric.egm = egm;

    userdata.electric.electrodeNames_uni = electrodeNames_uni;
    userdata.electric.egmUniX = unipolarEgmsX;
    userdata.electric.egmUni = unipolarEgms;

    userdata.electric.egmRefNames = nameRef;
    userdata.electric.egmRef = ref;
    userdata.electric.ecgNames = nameEcg;
    userdata.electric.ecg = squeeze(ecg);

    userdata.electric.annotations.woi = pointExport_WOI;
    userdata.electric.annotations.referenceAnnot = pointExport_ReferenceAnnotation;
    userdata.electric.annotations.mapAnnot = pointExport_MapAnnotation;
    userdata.electric.voltages.bipolar = pointExport_Bipolar;
    userdata.electric.voltages.unipolar = pointExport_Unipolar;
    userdata.electric.impedances.time = pointExport_ImpedanceTime;
    userdata.electric.impedances.value = pointExport_ImpedanceValue;


    %if ~isempty(filename_force)
    if nPoints>0
        if ~isempty(filename_force)
            userdata.electric.force.force = force;
            userdata.electric.force.axialAngle = axialAngle;
            userdata.electric.force.lateralAngle = lateralAngle;
            userdata.electric.force.time_force = t_force;
            userdata.electric.force.time_axial = t_axialAngle;
            userdata.electric.force.time_lateral = t_lateralAngle;
        end
    end

    % Now get the normal to the surface
    userdata.notes{1} = [date() ': Created'];
    [userdata.surface.triRep, userdata.surface.isVertexAtRim, userdata.surface.act_bip, normals, userdata.surface.uni_imp_frc] = read_meshfile([ studyDir, filesep(), map.meshFile]);
    if ~isempty(userdata.surface.triRep)
        warning('off', 'FINDCLOSESTVERTEX:hasToTrim')
        [closestVertices, ~] = findclosestvertex(userdata.surface.triRep, userdata.electric.egmX, true);
        warning('on', 'FINDCLOSESTVERTEX:hasToTrim')

        % Now work out the surface projections
        userdata.electric.egmSurfX = userdata.surface.triRep.X(closestVertices,:);
        userdata.electric.barDirection = normals(closestVertices, :);
    else
        userdata.notes{end+1,1} = [date ': No .mesh file found'];
        userdata.electric.egmSurfX = [];
        userdata.electric.barDirection = [];
    end

    % Now store the CF and RF data if the files existed
    if ~isempty(iContactForceInRfFiles)
        userdata.rf.originaldata.force.time = contactForceInRfData(:,1);
        userdata.rf.originaldata.force.force = contactForceInRfData(:,2);
        userdata.rf.originaldata.force.axialangle = contactForceInRfData(:,3);
        userdata.rf.originaldata.force.lateralangle = contactForceInRfData(:,4);

        X = NaN(size(contactForceInRfData(:,end),1),3); % was c - with NaNs in it
        b = contactForceInRfData(:,end);
        if ~isfield(userdata.electric, 'egmSurfX')
            warning('IMPORTCARTO_MEM: there was no meshfile, for userdata.rf.originaldata.force.position the position was used rather than position on shell')
            egmX = userdata.electric.egmX;
        else
            egmX = userdata.electric.egmSurfX; % was a
        end
        X(find(~isnan(b)),:) = egmX(b(find(~isnan(b))),:); % this is the line that does the location parsing into our matrix
        userdata.rf.originaldata.force.position = X;

    end
    if ~isempty(iRfFiles)
        userdata.rf.originaldata.ablparams.time = rfData(:,1);
        userdata.rf.originaldata.ablparams.power = rfData(:,3);
        userdata.rf.originaldata.ablparams.impedance = rfData(:,4);
        userdata.rf.originaldata.ablparams.distaltemp = rfData(:,5);
    end
    
    % restore old format - to be removed in future release
    for i = 1:numel(userdata.electric.electrodeNames_uni)
        userdata.electric.electrodeNames_uni{i} = [userdata.electric.electrodeNames_uni{i} , '('];
    end
    userdata.electric.egmRefNames = nameRefFull;
    userdata.electric.ecgNames = nameEcgFull;
    
    



    % Encourage user to save the data
    if ~isempty(saveFileName_cli)
        save(saveFileName_cli, 'userdata');
        matFileFullPath = saveFileName_cli;
    else
        defaultName = [map.studyName '_' map.name];
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


    %last of all remove the zip folder
end



% catch err
%     if ~isempty(hWait) && ishandle(hWait)
%         delete(hWait)
%     end
%     warning('IMPORTCARTO_MEM: Do you have extra Carto3 files on the path?');
%     rethrow(err);
% end

end

function fname = mycheckfilename(filename, allfilenames, searchstring)
% Check that filename has an exact match in allfilenames. If not, then
% search through filenames to see if there is a single string that contains
% searchstring. If not then return empty string;
tf = strcmp(filename, allfilenames);
if any(tf)
    fname = filename;
    return
end

k = strfind(allfilenames, searchstring);
iMatch = [];
for i = 1:length(k)
    if ~isempty(k{i})
        if ~isempty(iMatch)
            % we have already found a match

            warning('IMPORTCARTO3: there is more than one file that could represent this point.')
            fname = [];
            return;
        else
            iMatch = i;
        end
    end
end
if isempty(iMatch)
    fname = [];
    % use this for dbugging but can give excessive number of warnings
    % warning(['IMPORTCARTO3: the filename relating to ' char(39) searchstring char(39) ' is unexpected and no match was found.'])
else
    fname = allfilenames{iMatch};
    warning(['IMPORTCARTO3: the filename relating to ' char(39) searchstring char(39) ' is unexpected but a match was found - ' fname])
end
end


