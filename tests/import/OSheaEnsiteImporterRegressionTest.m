classdef OSheaEnsiteImporterRegressionTest < matlab.unittest.TestCase
    % Regression tests for malformed waveform exports in O'Shea Study2.

    properties
        TestDataRoot
    end

    methods (TestClassSetup)
        function locateTestData(testCase)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(repoRoot);

            testDataRoot = getenv('OPENEP_TESTING_DATA');
            if isempty(testDataRoot)
                candidate = fileparts(repoRoot);
                if isfolder(fullfile(candidate, 'EnsiteXDXLFiles'))
                    testDataRoot = candidate;
                end
            end
            testCase.TestDataRoot = testDataRoot;
        end
    end

    methods (Test)
        function waveRefsUnlabeledSignalColumnLoads(testCase)
            csvPath = fullfile(testCase.TestDataRoot, 'EnsiteXDXLFiles', ...
                'Study2-OShea-Birmingham', 'Wave_refs.csv');
            testCase.assumeTrue(isfile(csvPath), ...
                'O''Shea Wave_refs.csv test data not found.');

            warningState = warning('off', 'all');
            cleanupWarnings = onCleanup(@() warning(warningState)); %#ok<NASGU>
            [info, varnames, data] = loadensitex_dxldata( ...
                csvPath, 'ShowProgress', false);

            testCase.verifyEqual(info.sampleFreq, 2000);
            testCase.verifyEqual(size(data, 1), 1888);
            testCase.verifyEqual(varnames{end}, 'signals');
            testCase.verifyEqual(numel(data{1, end}), 1);
        end

        function waveRovMissingValuesPreserveSignalAlignment(testCase)
            csvPath = fullfile(testCase.TestDataRoot, 'EnsiteXDXLFiles', ...
                'Study2-OShea-Birmingham', 'Wave_rov.csv');
            testCase.assumeTrue(isfile(csvPath), ...
                'O''Shea Wave_rov.csv test data not found.');

            warningState = warning('off', 'all');
            cleanupWarnings = onCleanup(@() warning(warningState)); %#ok<NASGU>
            [info, varnames, data] = loadensitex_dxldata( ...
                csvPath, 'ShowProgress', false);

            testCase.verifyEqual(info.sampleFreq, 2000);
            testCase.verifyEqual(size(data, 1), 5991);
            testCase.verifyEqual(varnames{end}, 'signals');
            testCase.verifyEqual(numel(data{1, end}), 2001);
        end
    end
end
