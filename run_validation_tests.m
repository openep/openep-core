%% Run validation and importer regression tests

clear; clc;

repoRoot = fileparts(mfilename('fullpath'));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'validation'));

testFolders = {
    fullfile(repoRoot, 'tests', 'validation')
    fullfile(repoRoot, 'tests', 'import')
};

results = runtests(testFolders, 'IncludeSubfolders', true);
disp(table(results))
assertSuccess(results)
