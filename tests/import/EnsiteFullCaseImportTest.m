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
        function importsSelectedEgmTypeAndValidatesUserdata(testCase)
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
            egmType = getenv('OPENEP_FULL_ENSITEX_EGMTYPE');
            if isempty(egmType)
                egmType = 'bi';
            end

            outputFile = [tempname, '.mat'];
            cleanupObj = onCleanup(@() deleteIfPresent(outputFile));

            [userdata, savedFile] = importensitex_openep( ...
                caseRoot, ...
                'maptoread', mapName, ...
                'egmtype', egmType, ...
                'maptype', 'asegm', ...
                'showprogress', false, ...
                'savefilename', outputFile);

            testCase.verifyTrue(isfile(savedFile));
            report = validate_mapping_input(userdata, 'openep_userdata');
            testCase.verifyEqual(report.numFail, 0, report.summary);
            testCase.verifyTrue(any(strcmp({report.checks.id}, ...
                'openep.userdata.electric.present')));
        end
    end
end

function tf = runFullImporterTests()
tf = any(strcmpi(getenv('RUN_FULL_IMPORTER_SMOKE_TESTS'), ...
    {'1', 'true', 'yes'}));
end

function deleteIfPresent(filePath)
if isfile(filePath)
    delete(filePath);
end
end
