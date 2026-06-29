function results = import_ensitex_study3(caseRoot, outputRoot, egmTypes)
%IMPORT_ENSITEX_STUDY3 Import the three Study3 EnSiteX EGM configurations.
%
% results = import_ensitex_study3(caseRoot, outputRoot)
% results = import_ensitex_study3(caseRoot, outputRoot, {'bi', 'uni'})
%
% One OpenEP MAT file is created per requested EGM configuration.

repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'validation'));

if nargin < 1 || isempty(caseRoot)
    caseRoot = getenv('OPENEP_FULL_ENSITEX_CASE');
    if isempty(caseRoot)
        caseRoot = fullfile(fileparts(repoRoot), 'full_cases', 'EnsiteX', ...
            'Brussels', 'Study3-Gharaviri-Brussels');
    end
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile(pwd, 'openep_imports');
end
if nargin < 3 || isempty(egmTypes)
    egmTypes = {'bi', 'omni', 'uni'};
end

caseRoot = char(caseRoot);
outputRoot = char(outputRoot);
egmTypes = cellstr(egmTypes);
mapName = 'VoXel SR 1 ENDO';

assert(isfolder(caseRoot), 'EnSiteX Study3 folder not found: %s', caseRoot);
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

preflight = validate_mapping_input(caseRoot, 'ensitex_openep', ...
    'validationLevel', 'quick', ...
    'maxWaveFiles', 6);
fprintf('Pre-import validation: %s\n', preflight.summary);
assert(preflight.numFail == 0, 'Study3 pre-import validation failed.');

results = repmat(struct('egmType', '', 'outputFile', '', ...
    'elapsedSeconds', NaN, 'validationReport', struct()), numel(egmTypes), 1);

for i = 1:numel(egmTypes)
    egmType = lower(char(egmTypes{i}));
    outputFile = fullfile(outputRoot, ...
        sprintf('Study3_Gharaviri_Brussels_%s.mat', egmType));

    fprintf('\nImporting Study3 EGM type: %s\n', egmType);
    timer = tic;
    [userdata, savedFile] = importensitex_openep( ...
        caseRoot, ...
        'maptoread', mapName, ...
        'egmtype', egmType, ...
        'maptype', 'asegm', ...
        'showprogress', false, ...
        'savefilename', outputFile);
    elapsedSeconds = toc(timer);

    report = validate_mapping_input(userdata, 'openep_userdata');
    fprintf('%s import completed in %.1f seconds: %s\n', ...
        egmType, elapsedSeconds, report.summary);
    assert(report.numFail == 0, ...
        'Imported %s userdata failed validation.', egmType);

    results(i).egmType = egmType;
    results(i).outputFile = savedFile;
    results(i).elapsedSeconds = elapsedSeconds;
    results(i).validationReport = report;
end
end
