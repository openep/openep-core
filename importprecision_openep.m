function [userdata, matFileFullPath] = importprecision_openep(varargin)
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
% IMPORTCARTO_MEM accepts the following parameter-value pairs
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
%
% Example of command line entry ...
%       userdata = importprecision_openep(<path to folder>, ...
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
    dialog_title = 'Select the folder containing the exported Precision files';
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
    error('IMPORTPRECISION_OPENEP: Only specify a folder as user input')
else
    studyDir = userinput;
    homeDir = studyDir;
end

%% Parse command line input
nStandardArgs = 1; % UPDATE VALUE
channelRef_cli = '';
channelECG_cli = '';
saveFileName_cli = '';
if nargin > nStandardArgs
    for i = nStandardArgs+1:2:nargin
        switch varargin{i}
            case 'refchannel'
                channelRef_cli = varargin{i+1};
            case 'ecgchannel'
                channelECG_cli = varargin{i+1};
                if ischar(channelECG_cli); channelECG_cli = {channelECG_cli}; end
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
geometryDir = latMapDir; % this is the same, by definition, as the latMapDir

% And including the set of ALL the automap folders
allSubFolds = nameFolds(studyDir);
automapDirs = allSubFolds(strstartcmpi('OpenEP_AutoMap', allSubFolds));

%% Parse the geometry and surface mapping data
% By loading the LAT XML file to get the geometry
% By loading the LAT XML file to get the LAT map (I know, unnecessary
% duplication but this gives us future flexibility)
% By loading the bipolar XML file to get the bipolar map
% By loading the unipolar XML file to get the unipolar map
data_geometry = loadprecision_modelgroups(fullfile(latMapDir{:}, 'DxLandmarkGeo.xml'));
TRI = data_geometry.dxgeo.triangles;
X = data_geometry.dxgeo.vertices(:,1);
Y = data_geometry.dxgeo.vertices(:,2);
Z = data_geometry.dxgeo.vertices(:,3);
tr = TriRep(TRI, X, Y, Z);
t.X = tr.X;
t.Triangulation = tr.Triangulation;
normals = data_geometry.dxgeo.normals;

data_latMap = loadprecision_modelgroups(fullfile(latMapDir{:}, 'DxLandmarkGeo.xml'));
data_bipolarVoltageMap = loadprecision_modelgroups(fullfile(bipMapDir{:}, 'DxLandmarkGeo.xml'));
data_unipolarVoltageMap = loadprecision_modelgroups(fullfile(uniMapDir{:}, 'DxLandmarkGeo.xml'));
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
dxldataBip = importprecision_dxldata(bipMapDir);
dxldataUni = importprecision_dxldata(uniMapDir);

%% Parse electrogram data
% By loading the automaps (not needed for Precision but might be needed for
% EnSiteX)

%% Save all files in the OpenEP structure

% General data
userdata = openep_createuserdata();
userdata.systemName = 'precision';
userdata.notes{1} = [date() ': Created'];
userdata.precisionFolder = studyDir;
userdata.electric.sampleFrequency = dxldataBip.sampleFreq;

% Geometry
userdata.surface.triRep = t;
surfaceData = data_geometry.dxgeo.surface_of_origin;
userdata = setSurfaceProperty(userdata, 'name', 'surfaceOfOrigin', 'map', surfaceData, 'definedOn', 'elements');
userdata.surface.normals = normals;

% Surface maps, removing invalid data beyond interpolation distance
userdata.surface.act_bip = [act bip];
userdata.surface.uni_imp_frc = [uni imp frc];

% Electric data
userdata.electric.electrodeNames_bip = dxldataBip.rovtrace_pts';
userdata.electric.egmX = [dxldataBip.rovingx', dxldataBip.rovingy', dxldataBip.rovingz'];
userdata.electric.egmSurfX = [dxldataBip.surfPtx', dxldataBip.surfPty', dxldataBip.surfPtz'];
userdata.electric.egmRef = dxldataBip.reftrace'; 
userdata.electric.egm = dxldataBip.rovtrace';
userdata.electric.annotations.referenceAnnot = dxldataBip.refLAT';
userdata.electric.annotations.mapAnnot = dxldataBip.rovLAT';
userdata.electric.annotations.woi = 1 - userdata.electric.annotations.referenceAnnot;
userdata.electric.annotations.woi(:,2) = size(userdata.electric.egm,2) - userdata.electric.annotations.referenceAnnot;
userdata.electric.voltages.bipolar = dxldataBip.peak2peak';
userdata.electric.include = dxldataBip.utilized';
userdata.electric.names = strcat('P', strsplit(num2str(dxldataBip.ptnumber)))';

userdata.electric.electrodeNames_uni = dxldataUni.rovtrace_pts';
userdata.electric.egmX = [dxldataUni.rovingx', dxldataUni.rovingy', dxldataUni.rovingz'];
userdata.electric.egmUniSurfX = [dxldataUni.surfPtx', dxldataUni.surfPty', dxldataUni.surfPtz'];
userdata.electric.egmUni = dxldataUni.rovtrace';
userdata.electric.egmUni(:,:,2) = 0; % since we only get one unipole channel from Precision
userdata.electric.annotations.referenceAnnotUni = dxldataUni.refLAT';
userdata.electric.annotations.mapAnnotUni = dxldataUni.rovLAT';
userdata.electric.voltages.unipolar = dxldataUni.peak2peak';

% Temp - remote signalMaps which, if empty, prevents the file being loaded
% in EP Workbench
userdata.surface = rmfield(userdata.surface, 'signalMaps');
userdata.electric.tags = cell(length(userdata.electric.names),1);

% Encourage user to save the data
matFileFullPath = [];
if ~isempty(saveFileName_cli)
    save(saveFileName_cli, 'userdata');
    matFileFullPath = saveFileName_cli;
else
    defaultName = [dxldataBip.study '_' dxldataBip.mapId];
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

end