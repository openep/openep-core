function [ userdata ] = import_opencarp(varargin)
% TEMPLATE Brief description goes here
%
% Usage:
%   [ vol ] = myFunction( in1, in2 )
% Where:
%   in1  - the input
%   in2  - the input
%   out  - the output
%
% MYFUNCTION Detailed description goes here
%
% Author: Steven Williams (2016)
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

persistent saveDir homeDir

matFileFullPath = [];

if isempty(saveDir) || ~ischar(saveDir)
    saveDir = 'C:\Users\Nick (07989 436 479)\NickData\Research\CartoExport\CartoExportProcessed';
end
if ~isdir(saveDir)
    if exist('stevematlabroot.m','file')==2
        saveDir = stevematlabroot();
    elseif exist('nickmatlabroot.m','file')==2
        saveDir = nickmatlabroot();
    else
        saveDir = matlabroot();
    end
end

if ~isdir(saveDir)
    saveDir = cd;
end

if isempty(homeDir)
    homeDir = 'C:\Users\Nick (07989 436 479)\NickData\Research\CartoExport';
end
if ~isdir(homeDir); homeDir = saveDir; end

userdata = [];
hWait = [];

if nargin >= 1
    userinput = varargin{1};
else
    dialog_title = 'Select the openCARP .igb export file';
    filterSpec = {'*.igb;*.mat;*.zip', 'Appropriate files (*.igb;*.mat;*.zip)' ; '*.igb','IGB (*.igb)' ; '*.zip','Zip (*.zip)' ; '*.mat','Matlab (*.mat)' ; '*.*','All files (*.*)'};
    if ispc()
        uiwait(msgbox(dialog_title,'modal'));
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

% try
studyDirInfo = dir(studyDir);
allfilenames = cell(length(studyDirInfo),1);
for i = 1:length(studyDirInfo)
    allfilenames{i} = studyDirInfo(i).name;
end

disp(['Accessing: ' userinput]);
hWait = waitbar(0, 'Loading data');
fid = fopen(userinput);
Volt = fread(fid, 'float');
Volt = Volt(275:end); % remove header lines
Voltstore = Volt;
delete(hWait);

% Load the electrogram location mesh files
hWait = waitbar(0, 'Looking for electrogram locations');
[~, meshfilename, ~] = fileparts(userinput);
ptsFilename = [meshfilename '.pts'];
ptsFilename = mycheckfilename(ptsFilename, allfilenames, studyDir);
if isempty(ptsFilename)
    error('OPENEP/IMPORT_OPENCARP: Unable to locate electrogram location .pts file.')
else
    waitbar(0.5);
    elemfilename = [meshfilename '.elem'];
    elemfilename = mycheckfilename(elemfilename, allfilenames, studyDir);
    if isempty(elemfilename)
        error('OPENEP/IMPORT_OPENCARP: Unable to locate electrogram location .elem file.');
    else
        waitbar(1);
    end
end
delete(hWait);

hWait = waitbar(0, 'Reading electrogram location files');
VerticesBilayer = dlmread(ptsFilename, ' ', 1, 0);
waitbar(0.5);
FacesBilayer = dlmread(elemfilename, ' ', 1, 1);
waitbar(1)
FacesBilayer = FacesBilayer(:, 1:3)+1;
FV_egm.X = VerticesBilayer/1000; % divide by 1000 as OpenEP set up to work in mm
FV_egm.Triangulation = FacesBilayer;
delete(hWait);

% Load the surface mesh files
hWait = waitbar(0, 'Looking for surface mesh files');
%surfaceMeshName = 'Labelled';
[~, surfaceMeshName, ~] = fileparts(userinput); % defined on electrograms
surfacePtsFilename = [surfaceMeshName '.pts'];
surfacePtsFilename = mycheckfilename(surfacePtsFilename, allfilenames, studyDir);
if isempty(surfacePtsFilename)
    error('OPENEP/IMPORT_OPENCARP: Unable to locate surface .pts file.')
else
    waitbar(0.5);
    surfaceElemFilename = [surfaceMeshName '.elem'];
    surfaceElemFilename = mycheckfilename(surfaceElemFilename, allfilenames, studyDir);
    if isempty(surfaceElemFilename)
        error('OPENEP/IMPORT_OPENCARP: Unable to locate mesh .elem file.');
    else
        waitbar(1);
    end
end
delete(hWait);

hWait = waitbar(0, 'Reading surface mesh files');
VerticesBilayer = dlmread(surfacePtsFilename, ' ', 1, 0);
waitbar(0.5);
FacesBilayer = dlmread(surfaceElemFilename, ' ', 1, 1);
waitbar(1)
FacesBilayer = FacesBilayer(:, 1:3)+1;
FV_surface.X = VerticesBilayer/1000; % divide by 1000 as OpenEP set up to work in mm
FV_surface.Triangulation = FacesBilayer;
delete(hWait);

hWait = waitbar(0, 'Formating electrogram voltage data');
NOPOINTS=length(FV_egm.X);
Volt=reshape(Voltstore(1:floor(numel(Voltstore)/NOPOINTS)*NOPOINTS),NOPOINTS,floor(numel(Voltstore)/NOPOINTS)); %space, time point
LAT=zeros(size(Volt, 1), 1);
voltSize = size(Volt,1);
for ind=1:voltSize
    DvDt=diff(Volt(ind, :));
    [~, loc]=min(DvDt);
    LAT(ind)=loc;
    waitbar(ind/voltSize);
end
delete(hWait);

% create userdata
userdata = openep_createuserdata();

% assign the surface
userdata = setMesh(userdata, FV_surface);

% assign the  voltage electrogram data
userdata.electric.egmX = FV_egm.X;
iV = findclosestvertex(FV_surface.X, FV_egm.X);
userdata.electric.egmSurfX = (FV_surface.X(iV,:));

userdata.electric.annotations.woi = [ones(size(Volt,1),1) repmat(size(Volt,2),size(Volt,1),1)];
userdata.electric.annotations.referenceAnnot = zeros(size(Volt,1),1);
userdata.electric.annotations.mapAnnot = LAT;

% set the sampling frequency
userdata.electric.sampleFrequency = 1;

disp('done');


    function fname = mycheckfilename(filename, allfilenames, studyDir)
        % check that filename has an exact match in allfilenames. If not,
        % check the directory above the location of filename. If not then
        % return an empty string;
        fname = [];
        tf = strcmp(filename, allfilenames);
        if any(tf)
            fname = fullfile(studyDir,filename);
            return
        end
        upStudyDirInfo = dir(fileparts(studyDir));
        allfilenames = cell(length(upStudyDirInfo),1);
        for iFile = 1:length(upStudyDirInfo)
            allfilenamesUp{iFile} = upStudyDirInfo(iFile).name;
        end
        tf = strcmp(filename,allfilenamesUp);
        if any(tf)
            fname = fullfile(fileparts(studyDir),filename);
            return
        end

    end


end