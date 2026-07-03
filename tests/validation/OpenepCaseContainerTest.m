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

        function invalidBipolarEgmLayoutFailsValidation(testCase)
            openepCase = createCase({'bi'});
            openepCase.datasets.userdata.electric.egmUni = zeros(2, 5, 3);

            report = validate_mapping_input(openepCase, 'openep_case');

            testCase.verifyTrue(any(strcmp({report.checks.id}, ...
                'openep.case.bi.egm_layout') & ...
                strcmp({report.checks.level}, 'fail')));
        end

        function invalidOmnipolarEgmLayoutFailsValidation(testCase)
            openepCase = createCase({'omni'});
            openepCase.datasets.userdata.electric.egmUniX = zeros(2, 3, 2);

            report = validate_mapping_input(openepCase, 'openep_case');

            testCase.verifyTrue(any(strcmp({report.checks.id}, ...
                'openep.case.omni.egm_layout') & ...
                strcmp({report.checks.level}, 'fail')));
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
    datasets(i).userdata = createMinimalUserdata(modes{i});
end

openepCase = struct();
openepCase.schemaName = 'OpenEP multi-dataset case';
openepCase.schemaVersion = '1.0';
openepCase.source = struct('system', 'ensitex', 'rootFolder', '/test');
openepCase.mapName = 'Test Map';
openepCase.datasets = datasets;
end

function userdata = createMinimalUserdata(mode)
userdata = struct();
userdata.surface = struct();
userdata.surface.triRep = struct();
userdata.surface.triRep.X = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
userdata.surface.triRep.Triangulation = [1 2 3; 1 3 4];
userdata.surface.act_bip = [0 1; 1 1; 2 1; 3 1];
userdata.electric = struct();
userdata.electric.egmX = [0.1 0.1 0; 0.8 0.1 0];
userdata.electric.egmSurfX = userdata.electric.egmX;
userdata.electric.egm = zeros(2, 5);
if strcmp(mode, 'bi')
    userdata.electric.egmUni = zeros(2, 5, 2);
    userdata.electric.egmUniX = zeros(2, 3, 2);
    userdata.electric.electrodeNames_uni = cell(2, 2);
elseif strcmp(mode, 'omni')
    userdata.electric.egmUni = zeros(2, 5, 3);
    userdata.electric.egmUniX = zeros(2, 3, 3);
    userdata.electric.electrodeNames_uni = cell(2, 3);
else
    userdata.electric.egmUni = userdata.electric.egm;
    userdata.electric.egmUniX = userdata.electric.egmX;
    userdata.electric.electrodeNames_uni = cell(2, 1);
end
userdata.electric.voltages = struct('bipolar', [1; 1]);
end

function deleteIfPresent(filePath)
if isfile(filePath)
    delete(filePath);
end
end
