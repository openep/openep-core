function manifest = discover_full_cases(caseRoot)
%DISCOVER_FULL_CASES Find local CARTO and EnSite full-case exports.
%
% manifest = discover_full_cases('./full_cases')
%
% The returned struct array is intentionally compact so tests and notebooks
% can choose representative cases without repeatedly scanning large folders.

if nargin < 1 || isempty(caseRoot)
    caseRoot = fullfile(pwd, 'full_cases');
end
caseRoot = char(caseRoot);

manifest = emptyManifest();
if ~isfolder(caseRoot)
    return
end

manifest = [manifest discoverCartoCases(caseRoot)]; %#ok<AGROW>
manifest = [manifest discoverEnsiteCases(caseRoot)]; %#ok<AGROW>
end

function cases = discoverCartoCases(caseRoot)
cases = emptyManifest();
cartoRoot = fullfile(caseRoot, 'Carto');
if ~isfolder(cartoRoot)
    return
end

zipFiles = visibleFiles(dir(fullfile(cartoRoot, '**', '*.zip')));
for i = 1:numel(zipFiles)
    zipPath = fullfile(zipFiles(i).folder, zipFiles(i).name);
    cases(end+1) = makeCase('carto', zipPath, true, zipPath, zipFiles(i).bytes, ... %#ok<AGROW>
        {}, {}, {}, 'CARTO ZIP export');
end

meshFiles = visibleFiles(dir(fullfile(cartoRoot, '**', '*.mesh')));
folders = unique({meshFiles.folder}, 'stable');
for i = 1:numel(folders)
    folderPath = folders{i};
    xmlFiles = visibleFiles(dir(fullfile(folderPath, '*.xml')));
    names = {xmlFiles.name};
    hasStudyXml = any(~contains(names, 'Point_Export') & ~contains(names, 'Points_Export'));
    if hasStudyXml
        cases(end+1) = makeCase('carto', folderPath, false, '', folderSizeBytes(folderPath), ... %#ok<AGROW>
            inferCartoMapNames(folderPath), {}, {}, 'Extracted CARTO export');
    end
end
end

function cases = discoverEnsiteCases(caseRoot)
cases = emptyManifest();
ensiteRoot = fullfile(caseRoot, 'EnsiteX');
if ~isfolder(ensiteRoot)
    return
end

modelFiles = visibleFiles(dir(fullfile(ensiteRoot, '**', 'Contact_Mapping_Model.xml')));
for i = 1:numel(modelFiles)
    exportFolder = modelFiles(i).folder;
    [candidateMaps, candidateMapFiles, egmTypes] = inferEnsiteMaps(exportFolder);
    cases(end+1) = makeCase('ensitex', exportFolder, false, '', ... %#ok<AGROW>
        folderSizeBytes(exportFolder), candidateMaps, candidateMapFiles, egmTypes, ...
        'Extracted EnSiteX export');
end
end

function names = inferCartoMapNames(folderPath)
meshFiles = visibleFiles(dir(fullfile(folderPath, '*.mesh')));
names = cell(1, numel(meshFiles));
for i = 1:numel(meshFiles)
    [~, names{i}] = fileparts(meshFiles(i).name);
end
names = unique(names, 'stable');
end

function [mapNames, mapFiles, egmTypes] = inferEnsiteMaps(exportFolder)
mapFilesInfo = visibleFiles(dir(fullfile(exportFolder, '**', 'Map_LAT_*.csv')));
mapNames = {};
mapFiles = {};
egmTypes = {};
for i = 1:numel(mapFilesInfo)
    filePath = fullfile(mapFilesInfo(i).folder, mapFilesInfo(i).name);
    mapName = readEnsiteHeaderToken(filePath, 'Map name\s*:\s*,\s*([^,\r\n]+)');
    if isempty(mapName)
        [~, mapName] = fileparts(filePath);
    end
    mapNames{end+1} = strrep(strtrim(mapName), sprintf('\t'), ' '); %#ok<AGROW>
    mapFiles{end+1} = filePath; %#ok<AGROW>

    [~, baseName] = fileparts(filePath);
    tokens = regexp(baseName, 'Map_LAT_(.+)$', 'tokens', 'once');
    if ~isempty(tokens)
        egmTypes{end+1} = tokens{1}; %#ok<AGROW>
    end
end
mapNames = unique(mapNames, 'stable');
egmTypes = unique(egmTypes, 'stable');
end

function value = readEnsiteHeaderToken(filePath, pattern)
value = '';
fid = fopen(filePath, 'r');
if fid == -1
    return
end
cleanupObj = onCleanup(@() fclose(fid));

for i = 1:120
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    tokens = regexp(line, pattern, 'tokens', 'once');
    if ~isempty(tokens)
        value = tokens{1};
        return
    end
end
end

function bytes = folderSizeBytes(folderPath)
files = dir(fullfile(folderPath, '**', '*'));
files = files(~[files.isdir]);
files = visibleFiles(files);
if isempty(files)
    bytes = 0;
else
    bytes = sum([files.bytes]);
end
end

function files = visibleFiles(files)
if isempty(files)
    return
end
names = {files.name};
files = files(~startsWith(names, '.') & ~startsWith(names, '._'));
end

function s = makeCase(caseType, path, isArchive, archivePath, sizeBytes, candidateMaps, candidateMapFiles, egmTypes, notes)
[~, name, ext] = fileparts(path);
s = struct();
s.caseType = caseType;
s.name = [name, ext];
s.path = path;
s.isArchive = isArchive;
s.archivePath = archivePath;
s.sizeBytes = sizeBytes;
s.candidateMaps = candidateMaps;
s.candidateMapFiles = candidateMapFiles;
s.egmTypes = egmTypes;
s.notes = notes;
end

function s = emptyManifest()
s = struct('caseType', {}, 'name', {}, 'path', {}, 'isArchive', {}, ...
    'archivePath', {}, 'sizeBytes', {}, 'candidateMaps', {}, ...
    'candidateMapFiles', {}, 'egmTypes', {}, 'notes', {});
end
