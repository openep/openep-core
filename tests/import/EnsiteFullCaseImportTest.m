classdef EnsiteFullCaseImportTest < matlab.unittest.TestCase
    % Opt-in integration test for a complete multi-folder EnSiteX case.

    methods (TestClassSetup)
        function addProjectPaths(~)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(repoRoot);
            addpath(fullfile(repoRoot, 'validation'));
        end
    end

    methods (Test)
        function importsAllModesIntoValidatedCaseContainer(testCase)
            testCase.assumeTrue(runFullImporterTests(), ...
                'Set RUN_FULL_IMPORTER_SMOKE_TESTS=1 to run full imports.');

            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            caseRoot = getenv('OPENEP_FULL_ENSITEX_CASE');
            if isempty(caseRoot)
                caseRoot = fullfile(fileparts(repoRoot), 'full_cases', ...
                    'EnsiteX', 'Brussels', 'Study3-Gharaviri-Brussels');
            end
            testCase.assumeTrue(isfolder(caseRoot), ...
                'Full EnSiteX Study3 test case was not found.');

            mapName = getenv('OPENEP_FULL_ENSITEX_MAP');
            if isempty(mapName)
                mapName = 'VoXel SR 1 ENDO';
            end
            outputFile = [tempname, '.mat'];
            cleanupObj = onCleanup(@() deleteConversionFiles(outputFile));

            result = convert_mapping_case(caseRoot, outputFile, ...
                'system', 'ensitex', ...
                'maptoread', mapName, ...
                'modes', {'bi', 'uni', 'omni'}, ...
                'validationlevel', 'standard');

            testCase.verifyTrue(result.success, result.error.message);
            testCase.verifyTrue(result.outputPublished);
            testCase.verifyTrue(isfile(outputFile));
            testCase.verifyTrue(isfile(result.statusFile));
            testCase.verifyTrue(isfile(result.logFile));
            testCase.verifyEqual(result.outputValidation.numFail, 0, ...
                result.outputValidation.summary);
            loaded = load(outputFile, 'openepCase');
            openepCase = loaded.openepCase;
            testCase.verifyEqual({openepCase.datasets.recordingMode}, ...
                {'bi', 'uni', 'omni'});

            expectedPoints = [4585, 3620, 7110];
            actualPoints = arrayfun(@(d) size(d.userdata.electric.egmX, 1), ...
                openepCase.datasets);
            testCase.verifyEqual(actualPoints(:)', expectedPoints);

            verifyEgmLayout(testCase, openepCase.datasets(1), 2);
            verifyEgmLayout(testCase, openepCase.datasets(2), 1);
            verifyEgmLayout(testCase, openepCase.datasets(3), 3);
            delete(cleanupObj);
        end
    end
end

function verifyEgmLayout(testCase, dataset, nComponents)
electric = dataset.userdata.electric;
nPoints = size(electric.egm, 1);
nSamples = size(electric.egm, 2);

if nComponents == 1
    testCase.verifySize(electric.egmUni, [nPoints, nSamples]);
    testCase.verifySize(electric.egmUniX, [nPoints, 3]);
else
    testCase.verifySize(electric.egmUni, ...
        [nPoints, nSamples, nComponents]);
    testCase.verifySize(electric.egmUniX, [nPoints, 3, nComponents]);
end
testCase.verifySize(electric.electrodeNames_uni, [nPoints, nComponents]);
end

function tf = runFullImporterTests()
tf = any(strcmpi(getenv('RUN_FULL_IMPORTER_SMOKE_TESTS'), ...
    {'1', 'true', 'yes'}));
end

function deleteConversionFiles(outputFile)
[folder, name] = fileparts(outputFile);
files = {
    outputFile
    fullfile(folder, [name, '.status.json'])
    fullfile(folder, [name, '.log.txt'])
    };
for i = 1:numel(files)
    if isfile(files{i})
        delete(files{i});
    end
end
end
