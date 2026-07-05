function result = convert_mapping_case(inputPath, outputFile, varargin)
%CONVERT_MAPPING_CASE Headless CARTO or EnSiteX conversion for integration.
%
% result = convert_mapping_case(inputPath, outputFile, ...
%     'system', 'carto', ...
%     'maptoread', '2-LA', ...
%     'refchannel', 'CS1-CS2', ...
%     'ecgchannel', 'V1');
%
% The MAT output is written only after input and OpenEP output validation.
% A JSON status file and text log are written on both success and failure.

p = inputParser;
addRequired(p, 'inputPath', @(x) ischar(x) || isstring(x));
addRequired(p, 'outputFile', @(x) ischar(x) || isstring(x));
addParameter(p, 'system', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'maptoread', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'modes', {}, @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(p, 'refchannel', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'ecgchannel', '', ...
    @(x) ischar(x) || isstring(x) || iscellstr(x));
addParameter(p, 'validationlevel', 'standard', ...
    @(x) ischar(x) || isstring(x));
addParameter(p, 'statusfilename', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'logfilename', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'throwonfailure', false, ...
    @(x) islogical(x) && isscalar(x));
parse(p, inputPath, outputFile, varargin{:});
opts = p.Results;

validationFolder = fullfile(fileparts(mfilename('fullpath')), 'validation');
if exist('validate_mapping_input', 'file') ~= 2
    addpath(validationFolder);
end

inputPath = char(inputPath);
outputFile = char(outputFile);
[outputFolder, outputName, outputExtension] = fileparts(outputFile);
if isempty(outputFolder)
    outputFolder = pwd;
    outputFile = fullfile(outputFolder, [outputName, outputExtension]);
end
if ~strcmpi(outputExtension, '.mat')
    error('convert_mapping_case:OutputExtension', ...
        'outputFile must have a .mat extension.');
end
ensureFolder(outputFolder);

statusFile = char(opts.statusfilename);
if isempty(statusFile)
    statusFile = fullfile(outputFolder, [outputName, '.status.json']);
end
logFile = char(opts.logfilename);
if isempty(logFile)
    logFile = fullfile(outputFolder, [outputName, '.log.txt']);
end
ensureParentFolder(statusFile);
ensureParentFolder(logFile);

result = emptyResult(inputPath, outputFile, statusFile, logFile);
totalStart = tic;
consoleLog = '';
failureException = [];
cleanupObj = onCleanup(@() []);

try
    sourceSystem = resolveSystem(inputPath, opts.system);
    result.sourceSystem = sourceSystem;

    preparationStart = tic;
    if strcmp(sourceSystem, 'carto')
        validateCartoSelections(opts);
        [preparedInput, cleanupObj, archiveInfo] = prepareCartoInput(inputPath);
        result.archive = archiveInfo;
    else
        preparedInput = inputPath;
    end
    result.timings.preparationSeconds = toc(preparationStart);

    validationStart = tic;
    result.inputValidation = validateInput( ...
        preparedInput, sourceSystem, opts);
    result.timings.inputValidationSeconds = toc(validationStart);
    if result.inputValidation.numFail > 0
        error('convert_mapping_case:InputValidationFailed', ...
            'Input validation failed: %s', result.inputValidation.summary);
    end

    lastwarn('');
    importStart = tic;
    [consoleLog, payload, importException] = evalc( ...
        'invokeImporter(preparedInput, sourceSystem, opts)');
    result.timings.importSeconds = toc(importStart);
    [warningMessage, warningId] = lastwarn();
    result.runtimeWarning = struct( ...
        'identifier', warningId, 'message', warningMessage);
    if ~isempty(importException)
        throw(importException);
    end

    validationStart = tic;
    result.outputValidation = validate_mapping_input( ...
        payload.value, payload.validationMode);
    result.timings.outputValidationSeconds = toc(validationStart);
    if result.outputValidation.numFail > 0
        error('convert_mapping_case:OutputValidationFailed', ...
            'OpenEP output validation failed: %s', ...
            result.outputValidation.summary);
    end

    saveStart = tic;
    savePayloadAtomically(payload, outputFile);
    result.timings.saveSeconds = toc(saveStart);
    result.outputPublished = true;
    result.success = true;
    if result.inputValidation.numWarning > 0 || ...
            result.outputValidation.numWarning > 0 || ...
            ~isempty(result.runtimeWarning.message)
        result.status = 'warning';
    else
        result.status = 'success';
    end
catch ME
    failureException = ME;
    result.success = false;
    result.status = 'failure';
    result.error = exceptionAsStruct(ME);
end

try
    delete(cleanupObj);
catch cleanupException
    result.cleanupWarning = struct( ...
        'identifier', cleanupException.identifier, ...
        'message', cleanupException.message);
    if result.success
        result.status = 'warning';
    end
end
result.finishedAt = timestampNow();
result.timings.totalSeconds = toc(totalStart);
writeTextAtomically(logFile, formatLog(result, consoleLog));
writeJsonAtomically(statusFile, result);

if ~result.success && opts.throwonfailure
    throw(failureException);
end
end

function result = emptyResult(inputPath, outputFile, statusFile, logFile)
result = struct();
result.schemaName = 'OpenEP mapping conversion result';
result.schemaVersion = '1.0';
result.success = false;
result.status = 'failure';
result.sourceSystem = '';
result.inputPath = inputPath;
result.outputFile = outputFile;
result.outputPublished = false;
result.statusFile = statusFile;
result.logFile = logFile;
result.startedAt = timestampNow();
result.finishedAt = '';
result.archive = struct();
result.inputValidation = struct();
result.outputValidation = struct();
result.runtimeWarning = struct('identifier', '', 'message', '');
result.cleanupWarning = struct('identifier', '', 'message', '');
result.error = struct('identifier', '', 'message', '', 'stack', []);
result.timings = struct( ...
    'preparationSeconds', 0, ...
    'inputValidationSeconds', 0, ...
    'importSeconds', 0, ...
    'outputValidationSeconds', 0, ...
    'saveSeconds', 0, ...
    'totalSeconds', 0);
end

function sourceSystem = resolveSystem(inputPath, requestedSystem)
sourceSystem = lower(strtrim(char(requestedSystem)));
if strcmp(sourceSystem, 'ensite')
    sourceSystem = 'ensitex';
end
if ~strcmp(sourceSystem, 'auto')
    if ~any(strcmp(sourceSystem, {'carto', 'ensitex'}))
        error('convert_mapping_case:UnknownSystem', ...
            'system must be auto, carto or ensitex.');
    end
    return
end

if isfile(inputPath)
    [~, ~, ext] = fileparts(inputPath);
    if any(strcmpi(ext, {'.zip', '.xml'}))
        sourceSystem = 'carto';
        return
    end
elseif isfolder(inputPath)
    if ~isempty(dir(fullfile(inputPath, '**', 'Contact_Mapping_Model.xml')))
        sourceSystem = 'ensitex';
        return
    end
    if ~isempty(dir(fullfile(inputPath, '*.mesh')))
        sourceSystem = 'carto';
        return
    end
end

error('convert_mapping_case:SystemDetectionFailed', ...
    'Could not identify input as CARTO or EnSiteX: %s', inputPath);
end

function validateCartoSelections(opts)
if isempty(opts.maptoread) || isempty(opts.refchannel) || isempty(opts.ecgchannel)
    error('convert_mapping_case:CartoSelectionsRequired', ...
        ['Headless CARTO conversion requires maptoread, refchannel ', ...
        'and ecgchannel.']);
end
end

function [studyXml, cleanupObj, archiveInfo] = prepareCartoInput(inputPath)
[~, ~, ext] = fileparts(inputPath);
if isfile(inputPath) && strcmpi(ext, '.xml')
    studyXml = inputPath;
    cleanupObj = onCleanup(@() []);
    archiveInfo = struct();
    return
end

[caseFolder, cleanupObj, archiveInfo] = prepare_carto_case(inputPath);
studyXml = findCartoStudyXml(caseFolder);
end

function studyXml = findCartoStudyXml(caseFolder)
xmlFiles = dir(fullfile(caseFolder, '*.xml'));
names = {xmlFiles.name};
isStudy = ~startsWith(names, '.') & ...
    ~contains(names, 'Point_Export') & ...
    ~contains(names, 'Points_Export');
xmlFiles = xmlFiles(isStudy);
if numel(xmlFiles) ~= 1
    error('convert_mapping_case:CartoStudyXml', ...
        'Expected one CARTO study XML in %s, found %d.', ...
        caseFolder, numel(xmlFiles));
end
studyXml = fullfile(xmlFiles(1).folder, xmlFiles(1).name);
end

function report = validateInput(preparedInput, sourceSystem, opts)
if strcmp(sourceSystem, 'carto')
    report = validate_mapping_input(fileparts(preparedInput), ...
        'carto_openep', ...
        'mapToRead', opts.maptoread, ...
        'refChannel', opts.refchannel, ...
        'validationLevel', opts.validationlevel);
else
    report = validate_mapping_input(preparedInput, ...
        'ensitex_openep', ...
        'mapToRead', opts.maptoread, ...
        'validationLevel', opts.validationlevel);
end
end

function [payload, caughtException] = invokeImporter( ...
        preparedInput, sourceSystem, opts)
payload = struct('variableName', '', 'validationMode', '', 'value', struct());
caughtException = [];
try
    if strcmp(sourceSystem, 'carto')
        userdata = importcarto_mem(preparedInput, ...
            'maptoread', opts.maptoread, ...
            'refchannel', opts.refchannel, ...
            'ecgchannel', opts.ecgchannel, ...
            'verbose', false);
        payload.variableName = 'userdata';
        payload.validationMode = 'openep_userdata';
        payload.value = userdata;
    else
        openepCase = importensitex_case(preparedInput, ...
            'maptoread', opts.maptoread, ...
            'modes', opts.modes, ...
            'showprogress', false);
        payload.variableName = 'openepCase';
        payload.validationMode = 'openep_case';
        payload.value = openepCase;
    end
catch ME
    caughtException = ME;
end
end

function savePayloadAtomically(payload, outputFile)
outputFolder = fileparts(outputFile);
temporaryFile = [tempname(outputFolder), '.mat'];
cleanupObj = onCleanup(@() deleteIfPresent(temporaryFile));

if strcmp(payload.variableName, 'userdata')
    userdata = payload.value;
    save(temporaryFile, 'userdata', '-v7.3');
elseif strcmp(payload.variableName, 'openepCase')
    openepCase = payload.value;
    save(temporaryFile, 'openepCase', '-v7.3');
else
    error('convert_mapping_case:UnknownPayload', ...
        'Importer returned an unsupported payload.');
end

[moved, message] = movefile(temporaryFile, outputFile, 'f');
if ~moved
    error('convert_mapping_case:OutputMoveFailed', ...
        'Could not finalize MAT output: %s', message);
end
delete(cleanupObj);
end

function value = exceptionAsStruct(exception)
stack = struct('file', {}, 'name', {}, 'line', {});
for i = 1:numel(exception.stack)
    stack(i) = struct( ...
        'file', exception.stack(i).file, ...
        'name', exception.stack(i).name, ...
        'line', exception.stack(i).line);
end
value = struct( ...
    'identifier', exception.identifier, ...
    'message', exception.message, ...
    'stack', stack);
end

function text = formatLog(result, consoleLog)
lines = {
    sprintf('OpenEP conversion status: %s', upper(result.status))
    sprintf('System: %s', result.sourceSystem)
    sprintf('Input: %s', result.inputPath)
    sprintf('Output: %s', result.outputFile)
    sprintf('Started: %s', result.startedAt)
    sprintf('Finished: %s', result.finishedAt)
    sprintf('Duration: %.3f seconds', result.timings.totalSeconds)
    };
if ~isempty(fieldnames(result.inputValidation))
    lines{end+1} = ['Input validation: ', result.inputValidation.summary];
end
if ~isempty(fieldnames(result.outputValidation))
    lines{end+1} = ['Output validation: ', result.outputValidation.summary];
end
if ~result.success
    lines{end+1} = sprintf('Error [%s]: %s', ...
        result.error.identifier, result.error.message);
end
lines{end+1} = '';
lines{end+1} = 'Importer output:';
lines{end+1} = consoleLog;
text = strjoin(lines, newline);
end

function writeJsonAtomically(filePath, value)
text = jsonencode(value, 'PrettyPrint', true);
writeTextAtomically(filePath, text);
end

function writeTextAtomically(filePath, text)
folder = fileparts(filePath);
temporaryFile = tempname(folder);
cleanupObj = onCleanup(@() deleteIfPresent(temporaryFile));
fid = fopen(temporaryFile, 'w');
if fid == -1
    error('convert_mapping_case:StatusWriteFailed', ...
        'Could not create temporary output: %s', temporaryFile);
end
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);
delete(fileCleanup);

[moved, message] = movefile(temporaryFile, filePath, 'f');
if ~moved
    error('convert_mapping_case:StatusMoveFailed', ...
        'Could not finalize %s: %s', filePath, message);
end
delete(cleanupObj);
end

function ensureParentFolder(filePath)
folder = fileparts(filePath);
if ~isempty(folder)
    ensureFolder(folder);
end
end

function ensureFolder(folder)
if ~isfolder(folder)
    [created, message] = mkdir(folder);
    if ~created
        error('convert_mapping_case:CreateFolderFailed', ...
            'Could not create output folder %s: %s', folder, message);
    end
end
end

function deleteIfPresent(filePath)
if isfile(filePath)
    delete(filePath);
end
end

function value = timestampNow()
value = char(datetime('now', ...
    'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX'));
end
