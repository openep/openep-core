function [openepCase, outputFile] = import_ensitex_study3(caseRoot, outputRoot, egmTypes)
%IMPORT_ENSITEX_STUDY3 Import Study3 into one multi-dataset MAT file.
%
% openepCase = import_ensitex_study3(caseRoot, outputRoot)
% openepCase = import_ensitex_study3(caseRoot, outputRoot, {'bi', 'uni'})

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

assert(isfolder(caseRoot), 'EnSiteX Study3 folder not found: %s', caseRoot);
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

outputFile = fullfile(outputRoot, 'Study3_Gharaviri_Brussels_all.mat');
[openepCase, outputFile] = importensitex_case( ...
    caseRoot, ...
    'maptoread', 'VoXel SR 1 ENDO', ...
    'modes', egmTypes, ...
    'showprogress', false, ...
    'savefilename', outputFile);

report = validate_mapping_input(openepCase, 'openep_case');
fprintf('Study3 case import: %s\nSaved to: %s\n', report.summary, outputFile);
assert(report.numFail == 0, 'Imported Study3 case failed validation.');
end
