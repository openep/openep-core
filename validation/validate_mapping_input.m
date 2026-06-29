function report = validate_mapping_input(inputPath, workflowMode, varargin)
%VALIDATE_MAPPING_INPUT Lightweight validation for CARTO and EnSite inputs.
%
% report = validate_mapping_input(inputPath, workflowMode)
%
% V1 scope:
%   - Stage 1: file presence and workflow compatibility.
%   - Stage 2: EnSite DXL/CARTO header and required-column checks.
%   - Stage 3: simple cross-file consistency checks.
%   - Stage 4: numeric sanity checks for coordinates and scalar fields.
%   - Stage 5: OpenEP userdata structure checks after import.

p = inputParser;
addRequired(p, 'inputPath', @(x) ischar(x) || isstring(x) || isstruct(x));
addRequired(p, 'workflowMode', @(x) ischar(x) || isstring(x));
addParameter(p, 'mapToRead', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'refChannel', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'egmType', 'bi', @(x) ischar(x) || isstring(x));
addParameter(p, 'validationLevel', 'standard', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'maxWaveFiles', 6, ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1 && x == floor(x));
parse(p, inputPath, workflowMode, varargin{:});

workflowMode = lower(char(workflowMode));
opts = p.Results;
opts.validationLevel = lower(char(opts.validationLevel));
validLevels = {'quick', 'standard', 'full'};

checks = emptyCheck();
report = struct();
if isstruct(inputPath)
    report.inputPath = '<OpenEP userdata struct>';
else
    inputPath = char(inputPath);
    report.inputPath = inputPath;
end
report.workflowMode = workflowMode;
report.validationLevel = opts.validationLevel;
report.maxWaveFiles = opts.maxWaveFiles;

if ~any(strcmp(opts.validationLevel, validLevels))
    checks = addCheck(checks, 'fail', 1, 'validation.level.invalid', ...
        sprintf('Unknown validationLevel: %s. Use quick, standard, or full.', opts.validationLevel), report.inputPath);
    report.checks = checks;
    report = finalizeReport(report);
    return
elseif isstruct(inputPath) && ~strcmp(workflowMode, 'openep_userdata')
    checks = addCheck(checks, 'fail', 1, 'workflow.input_type_invalid', ...
        'Struct input is only supported with workflowMode=''openep_userdata''.', report.inputPath);
    report.checks = checks;
    report = finalizeReport(report);
    return
elseif ~isstruct(inputPath) && ~(isfolder(inputPath) || isfile(inputPath))
    checks = addCheck(checks, 'fail', 1, 'path.missing', ...
        sprintf('Input path does not exist: %s', inputPath), inputPath);
    report.checks = checks;
    report = finalizeReport(report);
    return
end

switch workflowMode
    case 'openep_userdata'
        checks = validateOpenepUserdata(inputPath, checks, report.inputPath);
    case 'openep_mat'
        checks = validateOpenepMat(inputPath, checks);
    case 'carto_openep'
        checks = validateCartoFolder(inputPath, checks, opts);
    case 'ensitex_openep'
        checks = validateEnsiteFolder(inputPath, checks, opts, true);
    case 'ensitex_dxl'
        checks = validateEnsiteFolder(inputPath, checks, opts, false);
    otherwise
        checks = addCheck(checks, 'fail', 1, 'workflow.unknown', ...
            sprintf('Unknown workflow mode: %s', workflowMode), inputPath);
end

report.checks = checks;
report = finalizeReport(report);
end

function checks = validateOpenepMat(inputPath, checks)
if ~isfile(inputPath)
    checks = addCheck(checks, 'fail', 1, 'openep_mat.not_file', ...
        'OpenEP MAT workflow expects a .mat file.', inputPath);
    return
end

[~, ~, ext] = fileparts(inputPath);
if ~strcmpi(ext, '.mat')
    checks = addCheck(checks, 'fail', 1, 'openep_mat.extension', ...
        'OpenEP MAT workflow expects a .mat file.', inputPath);
    return
end

try
    contents = whos('-file', inputPath);
    hasUserdata = any(strcmp({contents.name}, 'userdata'));
    if hasUserdata
        checks = addCheck(checks, 'pass', 1, 'openep_mat.userdata', ...
            'MAT file contains userdata.', inputPath);
        loaded = load(inputPath, 'userdata');
        checks = validateOpenepUserdata(loaded.userdata, checks, inputPath);
    else
        checks = addCheck(checks, 'fail', 1, 'openep_mat.no_userdata', ...
            'MAT file does not contain a variable named userdata.', inputPath);
    end
catch ME
    checks = addCheck(checks, 'fail', 1, 'openep_mat.unreadable', ...
        ['Could not inspect MAT file: ', ME.message], inputPath);
end
end

function checks = validateOpenepUserdata(userdata, checks, sourceLabel)
if ~isstruct(userdata)
    checks = addCheck(checks, 'fail', 5, 'openep.userdata.not_struct', ...
        'OpenEP userdata is expected to be a struct.', sourceLabel);
    return
end

checks = addCheck(checks, 'pass', 5, 'openep.userdata.struct', ...
    'OpenEP userdata is a struct.', sourceLabel);

if ~isfield(userdata, 'surface') || ~isstruct(userdata.surface)
    checks = addCheck(checks, 'fail', 5, 'openep.userdata.surface.missing', ...
        'userdata.surface is missing or is not a struct.', sourceLabel);
    return
else
    checks = addCheck(checks, 'pass', 5, 'openep.userdata.surface.present', ...
        'userdata.surface is present.', sourceLabel);
end

if ~isfield(userdata, 'electric') || ~isstruct(userdata.electric)
    checks = addCheck(checks, 'warning', 5, 'openep.userdata.electric.missing', ...
        'userdata.electric is missing or is not a struct; point-level checks will be limited.', sourceLabel);
else
    checks = addCheck(checks, 'pass', 5, 'openep.userdata.electric.present', ...
        'userdata.electric is present.', sourceLabel);
end

[vertices, faces, meshError] = extractOpenepMesh(userdata);
if isempty(meshError)
    checks = addCheck(checks, 'pass', 5, 'openep.mesh.readable', ...
        sprintf('Mesh was read: %d vertices, %d faces.', size(vertices, 1), size(faces, 1)), sourceLabel);
    checks = validateCoordinateMatrix(checks, vertices, 'openep.numeric.mesh_coordinates', ...
        'Mesh vertex coordinates', sourceLabel, true);
    checks = validateMeshFaces(checks, faces, size(vertices, 1), sourceLabel);
else
    checks = addCheck(checks, 'fail', 5, 'openep.mesh.unreadable', meshError, sourceLabel);
end

nPoints = NaN;
if exist('getNumPts', 'file') == 2
    try
        nPoints = getNumPts(userdata);
        checks = addCheck(checks, 'pass', 5, 'openep.userdata.num_points', ...
            sprintf('OpenEP reports %d mapping point(s).', nPoints), sourceLabel);
    catch ME
        checks = addCheck(checks, 'warning', 5, 'openep.userdata.num_points_unavailable', ...
            ['Could not call getNumPts: ', ME.message], sourceLabel);
    end
elseif isfield(userdata, 'electric') && isstruct(userdata.electric) && isfield(userdata.electric, 'egmX')
    nPoints = size(userdata.electric.egmX, 1);
    checks = addCheck(checks, 'info', 5, 'openep.userdata.num_points_fallback', ...
        sprintf('Estimated %d mapping point(s) from userdata.electric.egmX.', nPoints), sourceLabel);
else
    checks = addCheck(checks, 'info', 5, 'openep.userdata.num_points_unavailable', ...
        'Mapping point count could not be estimated.', sourceLabel);
end

if isfield(userdata, 'electric') && isstruct(userdata.electric)
    if isfield(userdata.electric, 'egmX')
        checks = validateCoordinateMatrix(checks, userdata.electric.egmX, ...
            'openep.numeric.egm_coordinates', 'Mapping point coordinates', sourceLabel, false);
        checks = validatePointCount(checks, userdata.electric.egmX, nPoints, ...
            'openep.userdata.egm_coordinates.rows', 'userdata.electric.egmX', sourceLabel);
    else
        checks = addCheck(checks, 'warning', 5, 'openep.userdata.egm_coordinates.missing', ...
            'userdata.electric.egmX is missing.', sourceLabel);
    end

    if isfield(userdata.electric, 'egmSurfX')
        checks = validateCoordinateMatrix(checks, userdata.electric.egmSurfX, ...
            'openep.numeric.surface_mapping_coordinates', 'Surface-projected mapping point coordinates', sourceLabel, false);
        checks = validatePointCount(checks, userdata.electric.egmSurfX, nPoints, ...
            'openep.userdata.egm_surface_coordinates.rows', 'userdata.electric.egmSurfX', sourceLabel);
    end

    if isfield(userdata.electric, 'voltages') && isstruct(userdata.electric.voltages) && ...
            isfield(userdata.electric.voltages, 'bipolar')
        bipolar = userdata.electric.voltages.bipolar;
        checks = validateNumericVector(checks, bipolar, 'openep.numeric.voltage.finite', ...
            'Point bipolar voltage', sourceLabel, false);
        checks = validateNonnegativeVector(checks, bipolar, 'openep.numeric.voltage.nonnegative', ...
            'Point bipolar voltage', sourceLabel);
        checks = validatePointCount(checks, bipolar, nPoints, ...
            'openep.userdata.voltage.rows', 'userdata.electric.voltages.bipolar', sourceLabel);
    end
end

if isfield(userdata.surface, 'act_bip')
    actBip = userdata.surface.act_bip;
    if ~isnumeric(actBip) || size(actBip, 2) < 2
        checks = addCheck(checks, 'warning', 5, 'openep.userdata.surface_act_bip.invalid', ...
            'userdata.surface.act_bip should be numeric with at least LAT and bipolar-voltage columns.', sourceLabel);
    else
        checks = validateNumericVector(checks, actBip(:, 1), 'openep.numeric.lat.finite', ...
            'Surface LAT', sourceLabel, false);
        checks = validateNumericVector(checks, actBip(:, 2), 'openep.numeric.surface_voltage.finite', ...
            'Surface bipolar voltage', sourceLabel, false);
        checks = validateNonnegativeVector(checks, actBip(:, 2), 'openep.numeric.surface_voltage.nonnegative', ...
            'Surface bipolar voltage', sourceLabel);
        if ~isempty(vertices)
            checks = validatePointCount(checks, actBip, size(vertices, 1), ...
                'openep.userdata.surface_act_bip.rows', 'userdata.surface.act_bip', sourceLabel);
        end
    end
else
    checks = addCheck(checks, 'warning', 5, 'openep.userdata.surface_act_bip.missing', ...
        'userdata.surface.act_bip is missing; LAT/voltage surface plotting may be limited.', sourceLabel);
end
end

function checks = validateCartoFolder(inputPath, checks, opts)
if isfile(inputPath)
    [~, ~, ext] = fileparts(inputPath);
    if strcmpi(ext, '.zip')
        fileInfo = dir(inputPath);
        checks = addCheck(checks, 'pass', 1, 'carto.archive.detected', ...
            sprintf('Found CARTO ZIP export (%.2f GB).', fileInfo.bytes / 1e9), inputPath);
        checks = addCheck(checks, 'warning', 1, 'carto.archive.requires_extraction', ...
            'CARTO ZIP export must be extracted or handled by batch conversion before content validation.', inputPath);
    else
        checks = addCheck(checks, 'fail', 1, 'carto.not_folder', ...
            'CARTO workflow expects an extracted folder or ZIP export.', inputPath);
    end
    return
elseif ~isfolder(inputPath)
    checks = addCheck(checks, 'fail', 1, 'carto.not_folder', ...
        'CARTO workflow expects a folder.', inputPath);
    return
end

xmlFiles = visibleFiles(dir(fullfile(inputPath, '*.xml')));
studyXml = xmlFiles(~contains({xmlFiles.name}, 'Point_Export'));
if isempty(studyXml)
    checks = addCheck(checks, 'fail', 1, 'carto.study_xml.missing', ...
        'No CARTO study XML found at the export root.', inputPath);
else
    checks = addCheck(checks, 'pass', 1, 'carto.study_xml.present', ...
        sprintf('Found CARTO study XML: %s', studyXml(1).name), ...
        fullfile(studyXml(1).folder, studyXml(1).name));
end

meshFiles = visibleFiles(dir(fullfile(inputPath, '*.mesh')));
if isempty(meshFiles)
    checks = addCheck(checks, 'fail', 1, 'carto.mesh.missing', ...
        'No CARTO mesh file found.', inputPath);
else
    checks = addCheck(checks, 'pass', 1, 'carto.mesh.present', ...
        sprintf('Found %d mesh file(s).', numel(meshFiles)), inputPath);
end

pointFiles = visibleFiles(dir(fullfile(inputPath, '*Point_Export.xml')));
if isempty(pointFiles)
    checks = addCheck(checks, 'fail', 1, 'carto.points.missing', ...
        'No CARTO point export XML files found.', inputPath);
else
    checks = addCheck(checks, 'pass', 1, 'carto.points.present', ...
        sprintf('Found %d point export XML file(s).', numel(pointFiles)), inputPath);
end

ecgFiles = visibleFiles(dir(fullfile(inputPath, '*ECG_Export*.txt')));
if isempty(ecgFiles)
    checks = addCheck(checks, 'warning', 1, 'carto.ecg.missing', ...
        'No ECG export files found. ECG/reference validation may be limited.', inputPath);
else
    checks = addCheck(checks, 'pass', 1, 'carto.ecg.present', ...
        sprintf('Found %d ECG export file(s).', numel(ecgFiles)), inputPath);
end

mapToRead = char(opts.mapToRead);
if ~isempty(mapToRead)
    mapFound = any(contains({meshFiles.name}, mapToRead));
    if ~mapFound && ~isempty(studyXml)
        studyText = safeFileRead(fullfile(studyXml(1).folder, studyXml(1).name));
        mapFound = contains(studyText, mapToRead);
    end

    if mapFound
        checks = addCheck(checks, 'pass', 2, 'carto.map.present', ...
            sprintf('Selected map was found: %s', mapToRead), inputPath);
    else
        checks = addCheck(checks, 'warning', 2, 'carto.map.not_found', ...
            sprintf('Selected map was not found by name: %s', mapToRead), inputPath);
    end
end

refChannel = char(opts.refChannel);
if ~isempty(refChannel) && ~isempty(ecgFiles)
    refFound = anyFileContains(ecgFiles, refChannel, 10);
    if refFound
        checks = addCheck(checks, 'pass', 2, 'carto.ref_channel.present', ...
            sprintf('Reference channel found in ECG exports: %s', refChannel), inputPath);
    else
        checks = addCheck(checks, 'warning', 2, 'carto.ref_channel.not_found', ...
            sprintf('Reference channel was not found in first ECG exports: %s', refChannel), inputPath);
    end
end
end

function checks = validateEnsiteFolder(inputPath, checks, opts, requireModelXml)
validationLevel = lower(char(opts.validationLevel));
countRows = ~strcmp(validationLevel, 'quick');
runNumericChecks = ~strcmp(validationLevel, 'quick');
unreadableLevel = 'fail';
if strcmp(validationLevel, 'quick')
    unreadableLevel = 'warning';
end

checks = addCheck(checks, 'info', 1, 'validation.level', ...
    sprintf('Validation level is %s.', validationLevel), inputPath);

if ~isfolder(inputPath)
    checks = addCheck(checks, 'fail', 1, 'ensite.not_folder', ...
        'EnSite workflow expects a folder.', inputPath);
    return
end

modelFiles = visibleFiles(dir(fullfile(inputPath, '**', 'Contact_Mapping_Model.xml')));
if requireModelXml
    if isempty(modelFiles)
        checks = addCheck(checks, 'fail', 1, 'ensite.model_xml.missing', ...
            'Contact_Mapping_Model.xml was not found. Full EnSite mesh import is not possible.', inputPath);
    else
        checks = addCheck(checks, 'pass', 1, 'ensite.model_xml.present', ...
            'Contact_Mapping_Model.xml was found.', fullfile(modelFiles(1).folder, modelFiles(1).name));
    end
else
    if isempty(modelFiles)
        checks = addCheck(checks, 'info', 1, 'ensite.model_xml.not_required', ...
            'DXL-only folder: mesh import is not expected without Contact_Mapping_Model.xml.', inputPath);
    else
        checks = addCheck(checks, 'info', 1, 'ensite.model_xml.present_in_dxl_mode', ...
            'Contact_Mapping_Model.xml is present; full EnSite workflow may also be possible.', ...
            fullfile(modelFiles(1).folder, modelFiles(1).name));
    end
end

csvFiles = visibleFiles(dir(fullfile(inputPath, '**', '*.csv')));
if isempty(csvFiles)
    checks = addCheck(checks, 'fail', 1, 'ensite.csv.missing', ...
        'No EnSite DXL CSV files found.', inputPath);
    return
end

latMapFiles = csvFiles(startsWith({csvFiles.name}, 'Map_LAT_', 'IgnoreCase', true));
voltageMapFiles = csvFiles(startsWith({csvFiles.name}, 'Map_PP_', 'IgnoreCase', true));
mapFiles = uniqueFiles([latMapFiles(:); voltageMapFiles(:)]);
if isempty(latMapFiles)
    checks = addCheck(checks, 'fail', 1, 'ensite.map_lat.missing', ...
        'No Map_LAT_*.csv file found.', inputPath);
else
    checks = addCheck(checks, 'pass', 1, 'ensite.map_lat.present', ...
        sprintf('Found %d LAT map CSV file(s).', numel(latMapFiles)), inputPath);
end
if ~isempty(voltageMapFiles)
    checks = addCheck(checks, 'info', 1, 'ensite.map_pp.present', ...
        sprintf('Found %d voltage map CSV file(s).', numel(voltageMapFiles)), inputPath);
end

allWaveFiles = csvFiles(startsWith({csvFiles.name}, 'Wave_', 'IgnoreCase', true));
waveRov = csvFiles(strcmpi({csvFiles.name}, 'Wave_rov.csv'));
waveRefs = csvFiles(strcmpi({csvFiles.name}, 'Wave_refs.csv'));
if isempty(waveRov)
    checks = addCheck(checks, 'warning', 1, 'ensite.wave_rov.missing', ...
        'Wave_rov.csv was not found. Waveform validation/import may be incomplete.', inputPath);
end
if isempty(waveRefs)
    checks = addCheck(checks, 'warning', 1, 'ensite.wave_refs.missing', ...
        'Wave_refs.csv was not found. Reference waveform validation/import may be incomplete.', inputPath);
end

egmType = lower(char(opts.egmType));
checks = validateEgmTypeFiles(csvFiles, checks, egmType, inputPath);

mapInfo = [];
if ~isempty(mapFiles)
    for i = 1:numel(mapFiles)
        filePath = fullfile(mapFiles(i).folder, mapFiles(i).name);
        info = inspectDxlCsv(filePath, countRows);
        checks = validateDxlHeader(info, checks, unreadableLevel);
        checks = validateMapColumns(info, checks, runNumericChecks);
        if isempty(mapInfo) && startsWith(info.name, 'Map_LAT_', 'IgnoreCase', true)
            mapInfo = info;
        end
    end
end

waveInfo = struct([]);
switch validationLevel
    case 'quick'
        waveFiles = selectCoreFirstWaveFiles(allWaveFiles, opts.maxWaveFiles);
        checks = addCheck(checks, 'info', 1, 'ensite.quick.wave_sample', ...
            sprintf(['Quick validation selected %d of %d Wave_*.csv file(s) ', ...
            'using core-first deterministic sampling.'], ...
            numel(waveFiles), numel(allWaveFiles)), inputPath);
    case 'full'
        waveFiles = allWaveFiles(:);
        checks = addCheck(checks, 'info', 1, 'ensite.full.wave_files', ...
            sprintf('Full validation selected all %d Wave_*.csv file(s).', numel(waveFiles)), inputPath);
    otherwise
        waveFiles = uniqueFiles([waveRefs(:); waveRov(:)]);
end

for i = 1:numel(waveFiles)
    filePath = fullfile(waveFiles(i).folder, waveFiles(i).name);
    info = inspectDxlCsv(filePath, countRows);
    checks = validateDxlHeader(info, checks, unreadableLevel);
    checks = validateWaveColumns(info, checks);
    waveInfo = [waveInfo info]; %#ok<AGROW>
end

if ~isempty(mapInfo) && ~isempty(waveInfo)
    checks = validateEnsiteConsistency(mapInfo, waveInfo, checks);
end
end

function checks = validateEgmTypeFiles(csvFiles, checks, egmType, inputPath)
names = {csvFiles.name};
switch egmType
    case 'bi'
        expected = {'wave_rov.csv', 'wave_refs.csv'};
    case 'omni'
        expected = {'wave_rov.csv', 'wave_refs.csv'};
    case 'uni'
        expected = {'wave_rov.csv', 'wave_refs.csv'};
    otherwise
        checks = addCheck(checks, 'warning', 1, 'ensite.egmtype.unknown', ...
            sprintf('Unknown egmType: %s', egmType), inputPath);
        return
end

for i = 1:numel(expected)
    if any(strcmpi(names, expected{i}))
        checks = addCheck(checks, 'pass', 1, ['ensite.egmtype.', expected{i}], ...
            sprintf('Required %s file is present for egmType=%s.', expected{i}, egmType), inputPath);
    else
        checks = addCheck(checks, 'warning', 1, ['ensite.egmtype.', expected{i}, '.missing'], ...
            sprintf('Expected %s for egmType=%s was not found.', expected{i}, egmType), inputPath);
    end
end
end

function checks = validateDxlHeader(info, checks, unreadableLevel)
if nargin < 3
    unreadableLevel = 'fail';
end
if ~isempty(info.error)
    checks = addCheck(checks, unreadableLevel, 2, 'ensite.csv.unreadable', info.error, info.file);
    return
end

if isempty(info.exportFileVersion)
    checks = addCheck(checks, 'fail', 2, 'ensite.header.version_missing', ...
        'Export File Version was not found.', info.file);
elseif any(strcmp(info.exportFileVersion, {'10.0R', '10', '11'})) || startsWith(info.exportFileVersion, '11')
    checks = addCheck(checks, 'pass', 2, 'ensite.header.version_supported', ...
        sprintf('Supported export file version: %s', info.exportFileVersion), info.file);
else
    checks = addCheck(checks, 'warning', 2, 'ensite.header.version_unknown', ...
        sprintf('Unexpected export file version: %s', info.exportFileVersion), info.file);
end

if isempty(info.dataElement)
    checks = addCheck(checks, 'fail', 2, 'ensite.header.data_element_missing', ...
        'Export Data Element was not found.', info.file);
elseif any(strcmpi(info.dataElement, {'DxL', 'DXLData'}))
    checks = addCheck(checks, 'pass', 2, 'ensite.header.data_element_valid', ...
        sprintf('Export Data Element is valid: %s', info.dataElement), info.file);
else
    checks = addCheck(checks, 'fail', 2, 'ensite.header.data_element_invalid', ...
        sprintf('Invalid Export Data Element: %s', info.dataElement), info.file);
end

if isnan(info.dataStartRow)
    checks = addCheck(checks, 'fail', 2, 'ensite.header.data_start_missing', ...
        'Data starts in row was not found.', info.file);
else
    checks = addCheck(checks, 'pass', 2, 'ensite.header.data_start_present', ...
        sprintf('Data starts in row %d.', info.dataStartRow), info.file);
end
end

function checks = validateMapColumns(info, checks, runNumericChecks)
if nargin < 3
    runNumericChecks = true;
end
if isempty(info.columns) || ~startsWith(info.name, 'Map_', 'IgnoreCase', true)
    return
end

required = {'(Point #)', 'Freeze Grp #', 'surface x', 'surface y', 'surface z', ...
    'roving x', 'roving y', 'roving z'};
if contains(lower(info.name), 'lat')
    required{end+1} = 'LAT';
end

checks = requireColumns(checks, info, required, 'ensite.map_columns');

if ~isnan(info.numMappingPts) && ~isnan(info.dataRows)
    if info.numMappingPts == info.dataRows
        checks = addCheck(checks, 'pass', 3, 'ensite.map_rows.match_header', ...
            sprintf('Map row count matches header: %d.', info.dataRows), info.file);
    else
        checks = addCheck(checks, 'warning', 3, 'ensite.map_rows.mismatch_header', ...
            sprintf('Map rows (%d) do not match header mapping points (%d).', ...
            info.dataRows, info.numMappingPts), info.file);
    end
end

if runNumericChecks
    checks = validateEnsiteMapNumeric(info, checks);
else
    checks = addCheck(checks, 'info', 4, 'ensite.numeric.skipped_quick', ...
        'Numeric map checks skipped in quick validation level.', info.file);
end
end

function checks = validateWaveColumns(info, checks)
if isempty(info.columns) || ~startsWith(info.name, 'Wave_', 'IgnoreCase', true)
    return
end

required = {'Trace', 'Freeze Grp #', '(Point #)', 'startTime (abs)', 'rovTime (wave samples)'};
checks = requireColumns(checks, info, required, 'ensite.wave_columns');

hasNumericSignalHeader = any(~isnan(str2double(info.columns)));
hasUnlabeledSignal = ~isempty(info.columns) && isempty(strtrim(info.columns{end}));
if hasNumericSignalHeader
    checks = addCheck(checks, 'pass', 2, 'ensite.wave_signal.numbered', ...
        'Wave file has numbered signal columns.', info.file);
elseif hasUnlabeledSignal
    checks = addCheck(checks, 'warning', 2, 'ensite.wave_signal.unlabeled', ...
        'Wave file has an unlabeled signal column; importer should store it as signals.', info.file);
else
    checks = addCheck(checks, 'warning', 2, 'ensite.wave_signal.missing', ...
        'No numeric or unlabeled signal column was detected.', info.file);
end

if isnan(info.sampleFreq)
    checks = addCheck(checks, 'fail', 2, 'ensite.wave_sample_rate.missing', ...
        'Sample rate was not found for wave file.', info.file);
elseif info.sampleFreq > 0
    checks = addCheck(checks, 'pass', 2, 'ensite.wave_sample_rate.valid', ...
        sprintf('Sample rate is %.3f Hz.', info.sampleFreq), info.file);
else
    checks = addCheck(checks, 'fail', 2, 'ensite.wave_sample_rate.invalid', ...
        sprintf('Sample rate is invalid: %.3f.', info.sampleFreq), info.file);
end
end

function checks = validateEnsiteConsistency(mapInfo, waveInfo, checks)
for i = 1:numel(waveInfo)
    info = waveInfo(i);
    if strcmpi(info.name, 'Wave_rov.csv') && ~isnan(mapInfo.numMappingPts) && ~isnan(info.dataRows)
        if abs(info.dataRows - mapInfo.numMappingPts) <= 1
            checks = addCheck(checks, 'pass', 3, 'ensite.wave_rov.rows_match_map', ...
                'Wave_rov row count is consistent with map points.', info.file);
        else
            checks = addCheck(checks, 'warning', 3, 'ensite.wave_rov.rows_mismatch_map', ...
                sprintf('Wave_rov rows (%d) differ from map points (%d).', ...
                info.dataRows, mapInfo.numMappingPts), info.file);
        end
    end
end

sampleFreqs = [waveInfo.sampleFreq];
sampleFreqs = sampleFreqs(~isnan(sampleFreqs));
if numel(unique(sampleFreqs)) <= 1 && ~isempty(sampleFreqs)
    checks = addCheck(checks, 'pass', 3, 'ensite.wave_sample_rate.consistent', ...
        'Wave files have consistent sample frequency.', mapInfo.file);
elseif numel(sampleFreqs) > 1
    checks = addCheck(checks, 'warning', 3, 'ensite.wave_sample_rate.inconsistent', ...
        'Wave files have inconsistent sample frequency.', mapInfo.file);
end
end

function checks = validateEnsiteMapNumeric(info, checks)
if isempty(info.columns) || ~startsWith(info.name, 'Map_', 'IgnoreCase', true)
    return
end

[surfaceX, surfaceFound] = readNumericColumns(info, {'surface x', 'surface y', 'surface z'});
if all(surfaceFound)
    checks = validateCoordinateColumns(checks, surfaceX, 'ensite.numeric.surface_coordinates', ...
        'Surface coordinates', info.file);
end

[rovingX, rovingFound] = readNumericColumns(info, {'roving x', 'roving y', 'roving z'});
if all(rovingFound)
    checks = validateCoordinateColumns(checks, rovingX, 'ensite.numeric.roving_coordinates', ...
        'Roving coordinates', info.file);
end

if hasColumn(info.columns, 'LAT')
    [lat, ~] = readNumericColumns(info, {'LAT'});
    checks = validateNumericVector(checks, lat, 'ensite.numeric.lat.finite', ...
        'LAT', info.file, false);
end

voltageColumns = {'P-P', 'peak2peak', 'pp_Valong', 'unipoleMaxPP'};
for i = 1:numel(voltageColumns)
    if hasColumn(info.columns, voltageColumns{i})
        [voltage, ~] = readNumericColumns(info, voltageColumns(i));
        checks = validateNumericVector(checks, voltage, 'ensite.numeric.voltage.finite', ...
            ['Voltage column ', voltageColumns{i}], info.file, false);
        checks = validateNonnegativeVector(checks, voltage, 'ensite.numeric.voltage.nonnegative', ...
            ['Voltage column ', voltageColumns{i}], info.file);
        break
    end
end
end

function checks = validateCoordinateColumns(checks, values, idPrefix, label, file)
if isempty(values)
    checks = addCheck(checks, 'warning', 4, [idPrefix, '.empty'], ...
        [label, ' could not be read.'], file);
    return
end

finiteRows = all(isfinite(values), 2);
if all(finiteRows)
    checks = addCheck(checks, 'pass', 4, [idPrefix, '.finite'], ...
        sprintf('%s are finite for %d row(s).', label, size(values, 1)), file);
elseif any(finiteRows)
    checks = addCheck(checks, 'warning', 4, [idPrefix, '.finite'], ...
        sprintf('%s contain non-finite values in %d of %d row(s).', ...
        label, sum(~finiteRows), size(values, 1)), file);
else
    checks = addCheck(checks, 'fail', 4, [idPrefix, '.finite'], ...
        [label, ' do not contain any fully finite coordinate rows.'], file);
end

if any(finiteRows)
    if all(abs(values(finiteRows, :)) < 1e-12, 'all')
        checks = addCheck(checks, 'warning', 4, [idPrefix, '.all_zero'], ...
            [label, ' are all zero; plotting may need roving coordinates or mesh data.'], file);
    else
        checks = addCheck(checks, 'pass', 4, [idPrefix, '.not_all_zero'], ...
            [label, ' are not all zero.'], file);
    end
end
end

function [values, found] = readNumericColumns(info, columnNames)
columnNames = cellstr(columnNames);
found = false(1, numel(columnNames));
indices = NaN(1, numel(columnNames));
for i = 1:numel(columnNames)
    indices(i) = columnIndex(info.columns, columnNames{i});
    found(i) = ~isnan(indices(i));
end

values = NaN(0, numel(columnNames));
if any(~found) || isnan(info.dataStartRow)
    return
end

fid = openTextFile(info.file);
if fid == -1
    return
end
cleanupObj = onCleanup(@() fclose(fid));

for i = 1:info.dataStartRow
    if ~ischar(fgetl(fid))
        return
    end
end

row = 0;
while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    if isempty(strtrim(line)) || strcmpi(strtrim(line), 'EOF')
        continue
    end
    parts = splitCsvLine(line);
    row = row + 1;
    values(row, :) = NaN;
    for j = 1:numel(columnNames)
        idx = indices(j);
        if idx <= numel(parts)
            values(row, j) = str2double(strtrim(parts{idx}));
        end
    end
end
end

function idx = columnIndex(columns, columnName)
idx = find(strcmpi(strtrim(columns), columnName), 1);
if isempty(idx)
    idx = NaN;
end
end

function checks = validateCoordinateMatrix(checks, values, idPrefix, label, file, required)
level = 'warning';
if required
    level = 'fail';
end

if ~isnumeric(values) || ~ismatrix(values) || size(values, 2) ~= 3
    checks = addCheck(checks, level, 4, [idPrefix, '.shape'], ...
        sprintf('%s should be a numeric N-by-3 array.', label), file);
    return
end

checks = addCheck(checks, 'pass', 4, [idPrefix, '.shape'], ...
    sprintf('%s have shape %d-by-3.', label, size(values, 1)), file);
checks = validateCoordinateColumns(checks, double(values), idPrefix, label, file);
end

function checks = validateNumericVector(checks, values, checkId, label, file, required)
level = 'warning';
if required
    level = 'fail';
end

if ~isnumeric(values) || isempty(values)
    checks = addCheck(checks, level, 4, checkId, ...
        sprintf('%s should be numeric and non-empty.', label), file);
    return
end

values = double(values(:));
finiteMask = isfinite(values);
if all(finiteMask)
    checks = addCheck(checks, 'pass', 4, checkId, ...
        sprintf('%s values are finite (%d value(s)).', label, numel(values)), file);
elseif any(finiteMask)
    checks = addCheck(checks, 'warning', 4, checkId, ...
        sprintf('%s has %d non-finite value(s) out of %d.', ...
        label, sum(~finiteMask), numel(values)), file);
else
    checks = addCheck(checks, level, 4, checkId, ...
        sprintf('%s does not contain any finite values.', label), file);
end
end

function checks = validateNonnegativeVector(checks, values, checkId, label, file)
if ~isnumeric(values) || isempty(values)
    checks = addCheck(checks, 'warning', 4, checkId, ...
        sprintf('%s cannot be checked for sign because it is not numeric or is empty.', label), file);
    return
end

values = double(values(:));
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    checks = addCheck(checks, 'warning', 4, checkId, ...
        sprintf('%s has no finite values for sign check.', label), file);
elseif any(finiteValues < 0)
    checks = addCheck(checks, 'warning', 4, checkId, ...
        sprintf('%s contains %d negative finite value(s).', label, sum(finiteValues < 0)), file);
else
    checks = addCheck(checks, 'pass', 4, checkId, ...
        sprintf('%s values are non-negative where finite.', label), file);
end
end

function checks = validatePointCount(checks, values, expectedRows, checkId, label, file)
if isnan(expectedRows) || ~isnumeric(values)
    return
end

actualRows = size(values, 1);
if actualRows == expectedRows
    checks = addCheck(checks, 'pass', 5, checkId, ...
        sprintf('%s row count matches expected count: %d.', label, expectedRows), file);
else
    checks = addCheck(checks, 'warning', 5, checkId, ...
        sprintf('%s has %d row(s), expected %d.', label, actualRows, expectedRows), file);
end
end

function checks = validateMeshFaces(checks, faces, nVertices, file)
if ~isnumeric(faces) || ~ismatrix(faces) || size(faces, 2) ~= 3 || isempty(faces)
    checks = addCheck(checks, 'fail', 5, 'openep.mesh.faces.shape', ...
        'Mesh faces should be a non-empty numeric N-by-3 array.', file);
    return
end

checks = addCheck(checks, 'pass', 5, 'openep.mesh.faces.shape', ...
    sprintf('Mesh faces have shape %d-by-3.', size(faces, 1)), file);

faces = double(faces);
validFaces = all(isfinite(faces), 2) & all(faces == round(faces), 2) & ...
    all(faces >= 1, 2) & all(faces <= nVertices, 2);
if all(validFaces)
    checks = addCheck(checks, 'pass', 5, 'openep.mesh.faces.indices', ...
        'Mesh face indices are finite integer vertex indices.', file);
elseif any(validFaces)
    checks = addCheck(checks, 'warning', 5, 'openep.mesh.faces.indices', ...
        sprintf('Mesh has %d invalid face row(s) out of %d.', ...
        sum(~validFaces), size(faces, 1)), file);
else
    checks = addCheck(checks, 'fail', 5, 'openep.mesh.faces.indices', ...
        'No mesh faces contain valid vertex indices.', file);
end
end

function [vertices, faces, meshError] = extractOpenepMesh(userdata)
vertices = [];
faces = [];
meshError = '';

if exist('getMesh', 'file') == 2
    try
        mesh = getMesh(userdata, 'type', 'struct');
        vertices = double(mesh.X);
        faces = double(mesh.Triangulation);
        return
    catch ME
        meshError = ['getMesh could not read userdata.surface.triRep: ', ME.message];
    end
end

if ~isfield(userdata, 'surface') || ~isstruct(userdata.surface) || ~isfield(userdata.surface, 'triRep')
    meshError = 'userdata.surface.triRep is missing.';
    return
end

triRep = userdata.surface.triRep;
try
    if isstruct(triRep) && isfield(triRep, 'X') && isfield(triRep, 'Triangulation')
        vertices = double(triRep.X);
        faces = double(triRep.Triangulation);
        meshError = '';
    elseif isa(triRep, 'triangulation')
        vertices = double(triRep.Points);
        faces = double(triRep.ConnectivityList);
        meshError = '';
    elseif isa(triRep, 'TriRep')
        vertices = double(triRep.X);
        faces = double(triRep.Triangulation);
        meshError = '';
    elseif isempty(meshError)
        meshError = 'userdata.surface.triRep is not a supported mesh representation.';
    end
catch ME
    meshError = ['Could not read userdata.surface.triRep: ', ME.message];
end
end

function filesOut = selectCoreFirstWaveFiles(filesIn, maxFiles)
filesIn = uniqueFiles(filesIn);
if numel(filesIn) <= maxFiles
    filesOut = filesIn;
    return
end

names = {filesIn.name};
isCore = strcmpi(names, 'Wave_refs.csv') | strcmpi(names, 'Wave_rov.csv');
coreFiles = filesIn(isCore);
otherFiles = filesIn(~isCore);

if numel(coreFiles) >= maxFiles
    filesOut = coreFiles(1:maxFiles);
    return
end

remaining = maxFiles - numel(coreFiles);
otherFiles = sampleFilesEvenly(otherFiles, remaining);
filesOut = uniqueFiles([coreFiles(:); otherFiles(:)]);
end

function filesOut = sampleFilesEvenly(filesIn, maxFiles)
filesIn = uniqueFiles(filesIn);
if isempty(filesIn) || maxFiles <= 0
    filesOut = filesIn([]);
    return
end
if numel(filesIn) <= maxFiles
    filesOut = filesIn;
    return
end

indices = unique(round(linspace(1, numel(filesIn), maxFiles)));
filesOut = filesIn(indices);
end

function filesOut = uniqueFiles(filesIn)
filesOut = filesIn;
if isempty(filesIn)
    return
end
paths = arrayfun(@(f) fullfile(f.folder, f.name), filesIn, 'UniformOutput', false);
[~, keep] = unique(paths, 'stable');
filesOut = filesIn(keep);
end

function checks = requireColumns(checks, info, requiredColumns, idPrefix)
missing = {};
for i = 1:numel(requiredColumns)
    if ~hasColumn(info.columns, requiredColumns{i})
        missing{end+1} = requiredColumns{i}; %#ok<AGROW>
    end
end

if isempty(missing)
    checks = addCheck(checks, 'pass', 2, [idPrefix, '.present'], ...
        'Required columns are present.', info.file);
else
    checks = addCheck(checks, 'fail', 2, [idPrefix, '.missing'], ...
        ['Missing required columns: ', strjoin(missing, ', ')], info.file);
end
end

function info = inspectDxlCsv(filePath, countRows)
if nargin < 2
    countRows = true;
end
[~, name, ext] = fileparts(filePath);
info = struct('file', filePath, 'name', [name, ext], 'error', '', ...
    'exportFileVersion', '', 'dataElement', '', 'mapName', '', 'mapType', '', ...
    'dataStartRow', NaN, 'numMappingPts', NaN, 'numFreezeGroups', NaN, ...
    'sampleFreq', NaN, 'waveformSamples', NaN, 'columns', {{}}, 'dataRows', NaN);

fid = openTextFile(filePath);
if fid == -1
    info.error = ['Could not open CSV file: ', filePath];
    return
end
cleanupObj = onCleanup(@() fclose(fid));

lines = {};
for i = 1:500
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    lines{end+1} = line; %#ok<AGROW>
    if startsWith(line, '*****') && ~isnan(info.dataStartRow) && numel(lines) >= info.dataStartRow
        break
    end
end

for i = 1:numel(lines)
    line = lines{i};
    info.exportFileVersion = firstToken(line, 'Export File Version\s*:\s*([^,\r\n]+)', info.exportFileVersion);
    info.dataElement = firstToken(line, 'Export Data Element\s*:\s*([^,\r\n]+)', info.dataElement);
    info.mapName = firstToken(line, 'Map name\s*:\s*,\s*([^,\r\n]+)', info.mapName);
    info.mapType = firstToken(line, 'Map type\s*:\s*,\s*([^,\r\n]+)', info.mapType);
    info.dataStartRow = firstNumber(line, 'Data starts in row\s*,\s*(\d+)', info.dataStartRow);
    info.numMappingPts = firstNumber(line, '# mapping pts\s*:\s*,\s*(\d+)', info.numMappingPts);
    info.numFreezeGroups = firstNumber(line, '# freeze groups\s*:\s*,\s*(\d+)', info.numFreezeGroups);
    info.sampleFreq = firstNumber(line, 'Sample rate\s*:\s*,\s*([0-9.]+)', info.sampleFreq);
    info.waveformSamples = firstNumber(line, 'Waveform samples exported\s*:\s*,\s*(\d+)', info.waveformSamples);
end

if ~isnan(info.dataStartRow) && numel(lines) >= info.dataStartRow
    info.columns = splitCsvLine(lines{info.dataStartRow});
elseif ~isnan(info.dataStartRow)
    info.columns = readSpecificLine(filePath, info.dataStartRow);
end

if countRows && ~isnan(info.dataStartRow)
    info.dataRows = countDataRows(filePath, info.dataStartRow);
end
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

function columns = readSpecificLine(filePath, lineNumber)
columns = {};
fid = openTextFile(filePath);
if fid == -1
    return
end
cleanupObj = onCleanup(@() fclose(fid));
line = '';
for i = 1:lineNumber
    line = fgetl(fid);
    if ~ischar(line)
        return
    end
end
columns = splitCsvLine(line);
end

function nRows = countDataRows(filePath, dataStartRow)
nRows = 0;
fid = openTextFile(filePath);
if fid == -1
    nRows = NaN;
    return
end
cleanupObj = onCleanup(@() fclose(fid));

for i = 1:dataStartRow
    if ~ischar(fgetl(fid))
        return
    end
end

while true
    line = fgetl(fid);
    if ~ischar(line)
        break
    end
    line = strtrim(line);
    if isempty(line) || strcmpi(line, 'EOF')
        continue
    end
    nRows = nRows + 1;
end
end

function fid = openTextFile(filePath)
fid = -1;
for attempt = 1:3
    fid = fopen(filePath, 'r');
    if fid ~= -1
        return
    end
    if attempt < 3
        pause(0.25 * attempt);
    end
end
end

function columns = splitCsvLine(line)
columns = regexp(line, ',', 'split');
end

function tf = hasColumn(columns, columnName)
tf = any(strcmpi(strtrim(columns), columnName));
end

function files = visibleFiles(files)
if isempty(files)
    return
end
names = {files.name};
files = files(~startsWith(names, '.') & ~startsWith(names, '._'));
end

function text = safeFileRead(filePath)
try
    text = fileread(filePath);
catch
    text = '';
end
end

function tf = anyFileContains(files, pattern, maxFiles)
tf = false;
for i = 1:min(numel(files), maxFiles)
    text = safeFileRead(fullfile(files(i).folder, files(i).name));
    if contains(text, pattern)
        tf = true;
        return
    end
end
end

function checks = emptyCheck()
checks = struct('level', {}, 'stage', {}, 'id', {}, 'message', {}, 'file', {});
end

function checks = addCheck(checks, level, stage, id, message, file)
checks(end+1) = struct( ...
    'level', char(level), ...
    'stage', stage, ...
    'id', char(id), ...
    'message', char(message), ...
    'file', char(file));
end

function report = finalizeReport(report)
levels = {report.checks.level};
if any(strcmp(levels, 'fail'))
    status = 'fail';
elseif any(strcmp(levels, 'warning'))
    status = 'warning';
else
    status = 'pass';
end

report.status = status;
report.numFail = sum(strcmp(levels, 'fail'));
report.numWarning = sum(strcmp(levels, 'warning'));
report.numPass = sum(strcmp(levels, 'pass'));
report.numInfo = sum(strcmp(levels, 'info'));
report.summary = sprintf('%s: %d fail, %d warning, %d pass, %d info', ...
    upper(status), report.numFail, report.numWarning, report.numPass, report.numInfo);
end
