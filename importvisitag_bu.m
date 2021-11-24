function [visitag, matFileFullPath] = importvisitag(varargin)
% IMPORTVISITAG provides a data structure from carto visitag files.
% Usage
%   visitag = importvisitag(dirName)
%   visitag = imporvisitag()
% Where:
%   dirName is the directory with all of the files corresponding to WiseTag
%   visitag is a single data structure
%
% IMPORTVISITAG detailed description goes here.
%
% visitag structure ...
%   .originaldata
%   .tag
%       .X
%       .surfX
%       .FTI
%   .grid
%       .X
%       .surfX
%       .FTI
%
% Author: Steven Williams (2014) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
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
if isempty(saveDir) || ~ischar(saveDir)
    if ispc
        saveDir = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
    else
        saveDir = getenv('HOME');
    end
end
if ~isdir(saveDir);
    if exist('stevematlabroot.m','file')==2
        saveDir = stevematlabroot();
    else
        saveDir = matlabroot();
    end
end
if isempty(homeDir)
    homeDir = saveDir;
end
if nargin >= 1 && ~strcmpi(varargin{1}, 'openfile')
    userinput = varargin{1};
else
    dialog_title = 'Select the WiseTag AblationSites.txt file';
    filterSpec = '*.txt';
    [filename,pathname] = uigetfile(filterSpec, dialog_title, homeDir);
    userinput = fullfile(pathname, filename);
    if filename == 0
        return
    end
end
visitag = [];
userdata = [];
if nargin>=2
    userdata = varargin{2};
end
studyDir = fileparts(userinput);
homeDir = studyDir;

% read the data files
try
    ablationSites = dlmread(userinput, '', 1,0); % AblationSites.txt
catch
    warning('No Visitag data found');
    if ~isempty(userdata)
        visitag = userdata;
    end
    return;
end
positionsData = dlmread([studyDir filesep() 'PositionsData.txt'], '', 1,0); % PositionsData.txt
contactForceData = dlmread([studyDir filesep() 'ContactForceData.txt'], '', 1,0); % ContactForceData.txt
sitesData = dlmread([studyDir filesep() 'Sites.txt'], '', 1,0);

% assemble data
[~,loc] = ismember(positionsData(:,5),contactForceData(:,1));
forceDataInPositions = contactForceData(loc,3);

%TODO ablationDataInPositions...

% parse the ablation sites getting the FTI
for i = 1:size(ablationSites,1)
    iStart = find(positionsData(:,1)==ablationSites(i,3)); %FirstPosTimeStamp
    iEnd = find(positionsData(:,1)==ablationSites(i,5));   % LastPosTimeStamp
    
    forces = forceDataInPositions(iStart:iEnd);
    FTI(i) = sum(forces/60);
end

FTI2 = sitesData(:,6) .* sitesData(:,7);

positions = sitesData(:,3:5);

visitag.calculatedFTI = FTI';
visitag.averageFTI = FTI2;
visitag.X = positions;

% now map the visitags to userdata if it has been passed in
if ~isempty(userdata)
    % Now work out the surface projections
    warning('off', 'FINDCLOSESTVERTEX:hasToTrim')
    [closestVertices, ~] = findclosestvertex(getMesh(userdata), visitag.X, true);
    warning('on', 'FINDCLOSESTVERTEX:hasToTrim')
    visitag.surfX = getMesh(userdata).X(closestVertices,:);

    userdata.visitag = visitag;
    visitag = userdata;
end
