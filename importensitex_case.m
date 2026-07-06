function [openepCase, matFileFullPath] = importensitex_case(studyDir, varargin)
%IMPORTENSITEX_CASE Import multiple EnSiteX recording modes into one MAT file.
%
% openepCase = importensitex_case(studyDir)
% openepCase = importensitex_case(studyDir, ...
%     'maptoread', mapName, ...
%     'modes', {'bi', 'uni', 'omni'}, ...
%     'savefilename', outputFile);
% Optional progresscallback is invoked as callback(stage, fraction, message).

p = inputParser;
addRequired(p, 'studyDir', @(x) ischar(x) || isstring(x));
addParameter(p, 'maptoread', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'modes', {}, @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(p, 'maptype', 'asegm', @(x) ischar(x) || isstring(x));
addParameter(p, 'savefilename', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'showprogress', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'progresscallback', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
parse(p, studyDir, varargin{:});
opts = p.Results;

studyDir = char(studyDir);
reportProgress(opts.progresscallback, 'discovering_exports', 0, ...
    'Discovering EnSiteX exports.');
manifest = inspectensitex_export(studyDir);
assert(~isempty(manifest.exports), ...
    'IMPORTENSITEX_CASE: No EnSiteX exports found in %s.', studyDir);

[mapName, mapExports] = selectMapExports(manifest.exports, opts.maptoread);
requestedModes = normalizeRequestedModes(opts.modes, mapExports);
datasets = repmat(emptyDataset(), numel(requestedModes), 1);
reportProgress(opts.progresscallback, 'discovered_exports', 0.05, ...
    sprintf('Found %d requested recording mode(s).', numel(requestedModes)));

for iMode = 1:numel(requestedModes)
    mode = requestedModes{iMode};
    export = selectModeExport(mapExports, mode);
    fprintf('\nImporting EnSiteX map "%s", mode %s\n', ...
        normalizeMapName(mapName), mode);
    reportProgress(opts.progresscallback, ['mode_', mode], ...
        (iMode - 1) / numel(requestedModes), ...
        sprintf('Importing EnSiteX mode %s (%d of %d).', ...
        mode, iMode, numel(requestedModes)));

    [userdata, ~] = importensitex_openep( ...
        export.folder, ...
        'maptoread', normalizeMapName(mapName), ...
        'egmtype', mode, ...
        'maptype', char(opts.maptype), ...
        'showprogress', opts.showprogress, ...
        'saveoutput', false);

    datasets(iMode).id = mode;
    datasets(iMode).recordingMode = mode;
    datasets(iMode).mapName = normalizeMapName(mapName);
    datasets(iMode).sourceFolder = export.folder;
    datasets(iMode).detection = struct( ...
        'confidence', export.confidence, ...
        'evidence', {export.evidence}, ...
        'warnings', {export.warnings});
    datasets(iMode).userdata = userdata;
    reportProgress(opts.progresscallback, ['mode_', mode], ...
        iMode / numel(requestedModes), ...
        sprintf('Finished EnSiteX mode %s (%d of %d).', ...
        mode, iMode, numel(requestedModes)));
end

openepCase = struct();
openepCase.schemaName = 'OpenEP multi-dataset case';
openepCase.schemaVersion = '1.0';
openepCase.source = struct( ...
    'system', 'ensitex', ...
    'rootFolder', studyDir, ...
    'importedAt', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss Z')));
openepCase.mapName = normalizeMapName(mapName);
openepCase.datasets = datasets;

matFileFullPath = char(opts.savefilename);
if ~isempty(matFileFullPath)
    save(matFileFullPath, 'openepCase', '-v7.3');
end

function reportProgress(callback, stage, fraction, message)
if ~isempty(callback)
    callback(stage, fraction, message);
end
end
end

function [mapName, mapExports] = selectMapExports(exports, requestedMap)
normalizedNames = cellfun(@normalizeMapName, {exports.mapName}, ...
    'UniformOutput', false);
uniqueNames = unique(normalizedNames, 'stable');

if isempty(requestedMap)
    if ~isscalar(uniqueNames)
        error(['IMPORTENSITEX_CASE: Multiple maps found. Specify maptoread. ', ...
            'Available maps: %s'], strjoin(uniqueNames, ', '));
    end
    selectedName = uniqueNames{1};
else
    requestedMap = normalizeMapName(requestedMap);
    matches = strcmpi(uniqueNames, requestedMap);
    if ~any(matches)
        matches = startsWith(uniqueNames, requestedMap, 'IgnoreCase', true);
    end
    if sum(matches) ~= 1
        error('IMPORTENSITEX_CASE: maptoread did not identify exactly one map.');
    end
    selectedName = uniqueNames{matches};
end

mapExports = exports(strcmp(normalizedNames, selectedName));
mapName = mapExports(1).mapName;
end

function modes = normalizeRequestedModes(requestedModes, exports)
if isempty(requestedModes)
    detected = {exports.recordingMode};
    canonical = {'bi', 'uni', 'omni'};
    modes = canonical(ismember(canonical, detected));
    if isempty(modes)
        error(['IMPORTENSITEX_CASE: No recording modes could be detected. ', ...
            'Specify modes explicitly.']);
    end
else
    modes = cellstr(requestedModes);
    modes = cellfun(@(x) lower(strtrim(x)), modes, 'UniformOutput', false);
    modes = unique(modes, 'stable');
end

validModes = {'bi', 'uni', 'omni'};
if any(~ismember(modes, validModes))
    error('IMPORTENSITEX_CASE: modes must contain only bi, uni or omni.');
end
end

function export = selectModeExport(exports, mode)
detectedModes = {exports.recordingMode};
matches = strcmp(detectedModes, mode);
if sum(matches) == 1
    export = exports(matches);
    return
end

unknown = strcmp(detectedModes, 'unknown');
if ~any(matches) && sum(unknown) == 1
    export = exports(unknown);
    return
end

error(['IMPORTENSITEX_CASE: Mode %s did not identify exactly one export. ', ...
    'Detected modes: %s'], mode, strjoin(detectedModes, ', '));
end

function name = normalizeMapName(name)
name = strrep(char(name), sprintf('\t'), ' ');
name = strtrim(regexprep(name, '\s+', ' '));
end

function dataset = emptyDataset()
dataset = struct( ...
    'id', '', ...
    'recordingMode', '', ...
    'mapName', '', ...
    'sourceFolder', '', ...
    'detection', struct(), ...
    'userdata', struct());
end
