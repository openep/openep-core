function manifest = inspectensitex_export(studyDir)
%INSPECTENSITEX_EXPORT Inspect EnSiteX exports without loading signal data.
%
% manifest = inspectensitex_export(studyDir)
%
% Folder and file names are treated as opaque. Semantic CSV headers are
% authoritative; recognized filenames are used only as a warned fallback.

studyDir = char(studyDir);
assert(isfolder(studyDir), 'EnSiteX study folder not found: %s', studyDir);

csvEntries = visibleFiles(dir(fullfile(studyDir, '**', '*.csv')));
files = emptyFileInfo();
ignoredFiles = {};
for i = 1:numel(csvEntries)
    filePath = fullfile(csvEntries(i).folder, csvEntries(i).name);
    info = inspectCsvHeader(filePath);
    if info.isDxl
        files(end+1) = info; %#ok<AGROW>
    else
        ignoredFiles{end+1} = filePath; %#ok<AGROW>
    end
end

manifest = struct();
manifest.studyDir = studyDir;
manifest.files = files;
manifest.ignoredFiles = ignoredFiles;
manifest.exports = groupExports(files);
end

function exports = groupExports(files)
exports = emptyExport();
if isempty(files)
    return
end

folders = unique({files.folder}, 'stable');
for iFolder = 1:numel(folders)
    inFolder = strcmp({files.folder}, folders{iFolder});
    folderFiles = files(inFolder);
    mapNames = unique({folderFiles.mapName}, 'stable');
    mapNames = mapNames(~cellfun('isempty', mapNames));
    for iMap = 1:numel(mapNames)
        inMap = strcmp({folderFiles.mapName}, mapNames{iMap});
        thisFiles = folderFiles(inMap);
        export = buildExport(thisFiles, numel(exports) + 1);
        exports(end+1) = export; %#ok<AGROW>
    end
end
end

function export = buildExport(files, index)
isMap = strcmp({files.kind}, 'map');
isWave = strcmp({files.kind}, 'wave');
mapFiles = files(isMap);
waveFiles = files(isWave);

[mode, confidence, evidence, warnings, errors] = detectMode(mapFiles, waveFiles);
numPoints = unique([mapFiles.numPoints]);
numPoints = numPoints(isfinite(numPoints));

geometryFile = findGeometryFile(files(1).folder);
export = struct();
export.id = sprintf('export_%d', index);
export.folder = files(1).folder;
export.mapName = files(1).mapName;
export.recordingMode = mode;
export.confidence = confidence;
export.evidence = evidence;
export.warnings = warnings;
export.errors = errors;
export.numPoints = numPoints;
export.geometryFile = geometryFile;
export.mapFiles = roleEntries(mapFiles);
export.waveFiles = roleEntries(waveFiles);
export.files = files;
end

function [mode, confidence, evidence, warnings, errors] = detectMode(mapFiles, waveFiles)
mode = 'unknown';
confidence = 'none';
evidence = {};
warnings = {};
errors = {};

tokens = {mapFiles.modeToken};
tokens = unique(tokens(~cellfun('isempty', tokens)));
if numel(tokens) > 1
    mode = 'conflict';
    errors{end+1} = ['Conflicting recording modes in Map type headers: ', ...
        strjoin(tokens, ', ')];
    return
elseif isscalar(tokens)
    mode = tokens{1};
    confidence = 'high';
    evidence{end+1} = ['Map type headers consistently identify ', mode, '.'];
    return
end

allColumns = {};
for i = 1:numel(mapFiles)
    allColumns = [allColumns mapFiles(i).columns]; %#ok<AGROW>
end
normalizedColumns = normalizeTokens(allColumns);
omniColumns = {'pp_vmax', 'pp_valong', 'pp_vacross', ...
    'uni_corner_elec', 'uni_along_elec', 'uni_across_elec'};
if all(ismember(omniColumns, normalizedColumns))
    mode = 'omni';
    confidence = 'high';
    evidence{end+1} = 'Omnipolar voltage and corner/along/across columns are present.';
    return
end

waveRoles = unique({waveFiles.role});
if all(ismember({'uni_corner', 'uni_along', 'uni_across'}, waveRoles))
    mode = 'omni';
    confidence = 'medium';
    evidence{end+1} = 'Corner, along and across unipolar wave roles are present.';
    return
end

filenameModes = unique([{mapFiles.filenameMode} {waveFiles.filenameMode}]);
filenameModes = filenameModes(~cellfun('isempty', filenameModes));
if isscalar(filenameModes)
    mode = filenameModes{1};
    confidence = 'low';
    evidence{end+1} = ['Legacy filenames identify ', mode, '.'];
    warnings{end+1} = 'Recording mode was inferred from filenames because semantic header evidence was absent.';
elseif numel(filenameModes) > 1
    errors{end+1} = ['Conflicting recording modes in legacy filenames: ', ...
        strjoin(filenameModes, ', ')];
end
end

function entries = roleEntries(files)
entryTemplate = struct('role', '', 'path', '', 'source', '');
entries = repmat(entryTemplate, numel(files), 1);
for i = 1:numel(files)
    entries(i) = struct( ...
        'role', files(i).role, ...
        'path', files(i).path, ...
        'source', files(i).roleSource);
end
end

function info = inspectCsvHeader(filePath)
[folder, name, ext] = fileparts(filePath);
info = emptyFileInfoScalar();
info.path = filePath;
info.folder = folder;
info.name = [name, ext];

fid = fopen(filePath, 'r');
if fid == -1
    info.error = 'Could not open file.';
    return
end
cleanupObj = onCleanup(@() fclose(fid));

lines = cell(1, 500);
nLines = 0;
for i = 1:500
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    nLines = nLines + 1;
    lines{nLines} = line;
    if nLines >= info.dataStartRow && isfinite(info.dataStartRow)
        break
    end
    info = parseHeaderLine(info, line);
end
lines = lines(1:nLines);

info.isDxl = any(strcmpi(info.dataElement, {'DxL', 'DXLData'}));
if ~info.isDxl
    return
end

if isfinite(info.dataStartRow) && numel(lines) >= info.dataStartRow
    info.columns = splitCsvLine(lines{info.dataStartRow});
end

if ~isempty(info.mapType) && ~strcmpi(info.mapType, 'N/A')
    info.kind = 'map';
    info.role = mapRole(info.mapType);
    info.roleSource = 'header';
    info.modeToken = modeFromToken(info.mapType);
else
    info.kind = 'wave';
    [info.role, info.roleSource] = waveRole(info.waveName, info.name);
end
info.filenameMode = modeFromFilename(info.name);
end

function info = parseHeaderLine(info, line)
info.dataElement = firstToken(line, ...
    'Export Data Element\s*:\s*([^,\r\n]+)', info.dataElement);
info.mapName = firstToken(line, ...
    'Map name\s*:\s*,\s*([^,\r\n]+)', info.mapName);
info.mapType = firstToken(line, ...
    'Map type\s*:\s*,\s*([^,\r\n]+)', info.mapType);
info.waveName = firstToken(line, ...
    'Wave name\s*:\s*,\s*([^,\r\n]+)', info.waveName);
info.dataStartRow = firstNumber(line, ...
    'Data starts in row\s*,\s*(\d+)', info.dataStartRow);
info.numPoints = firstNumber(line, ...
    '# mapping pts\s*:\s*,\s*(\d+)', info.numPoints);
info.numPoints = firstNumber(line, ...
    '# freeze groups\s*:\s*,\s*(\d+)', info.numPoints);
end

function role = mapRole(mapType)
token = regexprep(lower(strtrim(mapType)), '_(bi|uni|omni)$', '');
role = ['map_', normalizeToken(token)];
end

function [role, source] = waveRole(waveName, filename)
if ~isempty(waveName)
    role = normalizeToken(waveName);
    source = 'header';
    return
end

[~, baseName] = fileparts(filename);
baseName = regexprep(baseName, '^Wave_', '', 'ignorecase');
role = normalizeToken(baseName);
source = 'filename';
end

function mode = modeFromToken(value)
tokens = regexp(lower(strtrim(value)), '_(bi|uni|omni)$', 'tokens', 'once');
if isempty(tokens)
    mode = '';
else
    mode = tokens{1};
end
end

function mode = modeFromFilename(filename)
[~, baseName] = fileparts(filename);
mode = modeFromToken(baseName);
end

function geometryFile = findGeometryFile(folder)
geometryFile = '';
candidates = {
    fullfile(folder, 'Contact_Mapping_Model.xml')
    fullfile(fileparts(folder), 'Contact_Mapping_Model.xml')
};
for i = 1:numel(candidates)
    if isfile(candidates{i})
        geometryFile = candidates{i};
        return
    end
end
end

function token = normalizeToken(value)
token = lower(strtrim(char(value)));
token = regexprep(token, '[^a-z0-9]+', '_');
token = regexprep(token, '^_+|_+$', '');
end

function tokens = normalizeTokens(values)
tokens = cellfun(@normalizeToken, values, 'UniformOutput', false);
tokens = unique(tokens);
end

function value = firstToken(line, pattern, currentValue)
value = currentValue;
tokens = regexp(line, pattern, 'tokens', 'once');
if ~isempty(tokens)
    value = strtrim(tokens{1});
end
end

function value = firstNumber(line, pattern, currentValue)
value = currentValue;
tokens = regexp(line, pattern, 'tokens', 'once');
if ~isempty(tokens)
    value = str2double(tokens{1});
end
end

function columns = splitCsvLine(line)
columns = regexp(line, ',', 'split');
end

function files = visibleFiles(files)
if isempty(files)
    return
end
names = {files.name};
files = files(~startsWith(names, '.') & ~startsWith(names, '._'));
end

function files = emptyFileInfo()
files = repmat(emptyFileInfoScalar(), 0, 1);
end

function info = emptyFileInfoScalar()
info = struct( ...
    'path', '', ...
    'folder', '', ...
    'name', '', ...
    'isDxl', false, ...
    'kind', '', ...
    'dataElement', '', ...
    'mapName', '', ...
    'mapType', 'N/A', ...
    'waveName', '', ...
    'dataStartRow', NaN, ...
    'numPoints', NaN, ...
    'columns', {{}}, ...
    'modeToken', '', ...
    'filenameMode', '', ...
    'role', '', ...
    'roleSource', '', ...
    'error', '');
end

function exports = emptyExport()
exports = struct( ...
    'id', {}, ...
    'folder', {}, ...
    'mapName', {}, ...
    'recordingMode', {}, ...
    'confidence', {}, ...
    'evidence', {}, ...
    'warnings', {}, ...
    'errors', {}, ...
    'numPoints', {}, ...
    'geometryFile', {}, ...
    'mapFiles', {}, ...
    'waveFiles', {}, ...
    'files', {});
end
