classdef OpenepCaseContainerTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProjectPaths(~)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(repoRoot);
            addpath(fullfile(repoRoot, 'validation'));
        end
    end

    methods (Test)
        function validatesAndSelectsDatasets(testCase)
            openepCase = createCase({'bi', 'omni'});

            report = validate_mapping_input(openepCase, 'openep_case');
            [userdata, dataset] = select_openep_dataset(openepCase, 'omni');

            testCase.verifyEqual(report.numFail, 0, report.summary);
            testCase.verifyEqual(dataset.recordingMode, 'omni');
            testCase.verifyEqual(userdata.electric.egmX, ...
                openepCase.datasets(2).userdata.electric.egmX);
        end

        function duplicateModesFailValidation(testCase)
            openepCase = createCase({'bi', 'bi'});

            report = validate_mapping_input(openepCase, 'openep_case');

            testCase.verifyGreaterThan(report.numFail, 0);
            testCase.verifyTrue(any(strcmp({report.checks.id}, ...
                'openep.case.modes.invalid')));
        end

        function validatesSavedCaseMatFile(testCase)
            openepCase = createCase({'uni'});
            matFile = [tempname, '.mat'];
            cleanupObj = onCleanup(@() deleteIfPresent(matFile));
            save(matFile, 'openepCase', '-v7.3');

            report = validate_mapping_input(matFile, 'openep_mat');

            testCase.verifyEqual(report.numFail, 0, report.summary);
            testCase.verifyTrue(any(strcmp({report.checks.id}, ...
                'openep_mat.case')));
            delete(cleanupObj);
        end
    end
end

function openepCase = createCase(modes)
datasets = repmat(struct( ...
    'id', '', ...
    'recordingMode', '', ...
    'mapName', 'Test Map', ...
    'sourceFolder', '/test', ...
    'detection', struct(), ...
    'userdata', struct()), numel(modes), 1);
for i = 1:numel(modes)
    datasets(i).id = sprintf('%s_%d', modes{i}, i);
    datasets(i).recordingMode = modes{i};
    datasets(i).userdata = createMinimalUserdata();
end

openepCase = struct();
openepCase.schemaName = 'OpenEP multi-dataset case';
openepCase.schemaVersion = '1.0';
openepCase.source = struct('system', 'ensitex', 'rootFolder', '/test');
openepCase.mapName = 'Test Map';
openepCase.datasets = datasets;
end

function userdata = createMinimalUserdata()
userdata = struct();
userdata.surface = struct();
userdata.surface.triRep = struct();
userdata.surface.triRep.X = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
userdata.surface.triRep.Triangulation = [1 2 3; 1 3 4];
userdata.surface.act_bip = [0 1; 1 1; 2 1; 3 1];
userdata.electric = struct();
userdata.electric.egmX = [0.1 0.1 0; 0.8 0.1 0];
userdata.electric.egmSurfX = userdata.electric.egmX;
userdata.electric.voltages = struct('bipolar', [1; 1]);
end

function deleteIfPresent(filePath)
if isfile(filePath)
    delete(filePath);
end
end
