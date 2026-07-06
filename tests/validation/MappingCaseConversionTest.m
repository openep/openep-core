classdef MappingCaseConversionTest < matlab.unittest.TestCase
    % Tests for the headless integration-facing conversion contract.

    properties
        TempRoot
    end

    methods (TestClassSetup)
        function addProjectPaths(testCase)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(repoRoot);
            addpath(fullfile(repoRoot, 'validation'));
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
        end
    end

    methods (TestClassTeardown)
        function removeTemporaryFolder(testCase)
            if isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, 's');
            end
        end
    end

    methods (Test)
        function failureWritesMachineReadableArtifacts(testCase)
            outputFile = fullfile(testCase.TempRoot, 'failed.mat');

            result = convert_mapping_case( ...
                fullfile(testCase.TempRoot, 'missing.zip'), outputFile, ...
                'system', 'carto', ...
                'maptoread', '2-LA', ...
                'refchannel', 'CS1-CS2', ...
                'ecgchannel', 'V1');

            testCase.verifyFalse(result.success);
            testCase.verifyFalse(result.outputPublished);
            testCase.verifyEqual(result.status, 'failure');
            testCase.verifyFalse(isfile(outputFile));
            testCase.verifyTrue(isfile(result.statusFile));
            testCase.verifyTrue(isfile(result.logFile));
            testCase.verifyFalse(isfile(result.progressFile));

            savedStatus = jsondecode(fileread(result.statusFile));
            testCase.verifyFalse(savedStatus.success);
            testCase.verifyFalse(savedStatus.outputPublished);
            testCase.verifyEqual(savedStatus.status, 'failure');
            testCase.verifyEqual(savedStatus.error.identifier, ...
                'prepare_carto_case:MissingPath');
            testCase.verifySubstring(fileread(result.logFile), ...
                'OpenEP conversion status: FAILURE');
        end

        function throwOnFailureWritesArtifactsBeforeError(testCase)
            outputFile = fullfile(testCase.TempRoot, 'thrown.mat');
            statusFile = fullfile(testCase.TempRoot, 'thrown.json');
            logFile = fullfile(testCase.TempRoot, 'thrown.log');
            progressFile = fullfile(testCase.TempRoot, 'thrown.progress.json');

            call = @() convert_mapping_case( ...
                fullfile(testCase.TempRoot, 'missing.zip'), outputFile, ...
                'system', 'carto', ...
                'maptoread', '2-LA', ...
                'refchannel', 'CS1-CS2', ...
                'ecgchannel', 'V1', ...
                'statusfilename', statusFile, ...
                'logfilename', logFile, ...
                'progressfilename', progressFile, ...
                'throwonfailure', true);

            testCase.verifyError(call, 'prepare_carto_case:MissingPath');
            testCase.verifyTrue(isfile(statusFile));
            testCase.verifyTrue(isfile(logFile));
            testCase.verifyFalse(isfile(progressFile));
            testCase.verifyFalse(isfile(outputFile));
        end
    end
end
