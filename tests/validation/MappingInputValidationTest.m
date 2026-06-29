classdef MappingInputValidationTest < matlab.unittest.TestCase
    % Unit tests for mapping-input and imported-userdata validation.

    properties
        RepoRoot
        TempRoot
    end

    methods (TestClassSetup)
        function setupPathsAndTemporaryFolder(testCase)
            testCase.RepoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(testCase.RepoRoot);
            addpath(fullfile(testCase.RepoRoot, 'validation'));
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
        function openepUserdataStructurePasses(testCase)
            userdata = createMinimalOpenepUserdata();

            report = validate_mapping_input(userdata, 'openep_userdata');

            testCase.verifyEqual(report.numFail, 0, report.summary);
            testCase.verifyTrue(hasCheck(report, 'openep.mesh.readable'));
            testCase.verifyTrue(hasCheck(report, 'openep.numeric.lat.finite'));
            testCase.verifyTrue(hasCheck(report, 'openep.numeric.voltage.nonnegative'));
        end

        function missingSurfaceFails(testCase)
            userdata = struct('electric', struct());

            report = validate_mapping_input(userdata, 'openep_userdata');

            testCase.verifyEqual(report.status, 'fail', report.summary);
            testCase.verifyTrue(hasCheck(report, 'openep.userdata.surface.missing'));
        end

        function discoversCartoAndEnsiteCases(testCase)
            caseRoot = fullfile(testCase.TempRoot, 'full_cases');
            cartoFolder = fullfile(caseRoot, 'Carto', 'Study1');
            ensiteFolder = fullfile(caseRoot, 'EnsiteX', 'Study1', 'Export_bi');
            mkdir(cartoFolder);
            mkdir(fullfile(ensiteFolder, 'Contact_Mapping'));

            writeTextFile(fullfile(cartoFolder, 'study.xml'), '<study />');
            writeTextFile(fullfile(cartoFolder, 'map.mesh'), 'mesh');
            writeTextFile(fullfile(ensiteFolder, 'Contact_Mapping_Model.xml'), '<model />');
            writeTextFile(fullfile(ensiteFolder, 'Contact_Mapping', 'Map_LAT_bi.csv'), ...
                sprintf('Export File Version: 11\nMap name:,Test Map\n'));

            manifest = discover_full_cases(caseRoot);

            testCase.verifyEqual(numel(manifest), 2);
            testCase.verifyTrue(any(strcmp({manifest.caseType}, 'carto')));
            testCase.verifyTrue(any(strcmp({manifest.caseType}, 'ensitex')));
            ensiteCase = manifest(strcmp({manifest.caseType}, 'ensitex'));
            testCase.verifyEqual(ensiteCase.candidateMaps, {'Test Map'});
            testCase.verifyEqual(ensiteCase.egmTypes, {'bi'});
        end

        function preparesAndCleansCartoZip(testCase)
            sourceFolder = fullfile(testCase.TempRoot, 'carto_zip_source');
            mkdir(sourceFolder);
            writeTextFile(fullfile(sourceFolder, 'study.xml'), '<study />');
            writeTextFile(fullfile(sourceFolder, 'map.mesh'), 'mesh');
            zipFile = fullfile(testCase.TempRoot, 'carto_case.zip');
            zip(zipFile, sourceFolder);

            [caseFolder, cleanupObj, info] = prepare_carto_case_for_test(zipFile);

            testCase.verifyTrue(isfolder(caseFolder));
            testCase.verifyTrue(info.wasArchive);
            testCase.verifyTrue(isfile(fullfile(caseFolder, 'map.mesh')));
            delete(cleanupObj);
            testCase.verifyFalse(isfolder(info.extractionRoot));
        end
    end
end

function userdata = createMinimalOpenepUserdata()
userdata = struct();
userdata.surface = struct();
userdata.surface.triRep = struct();
userdata.surface.triRep.X = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
userdata.surface.triRep.Triangulation = [1 2 3; 1 3 4];
userdata.surface.act_bip = [0 1.1; 5 0.9; 10 1.4; 15 0.8];
userdata.electric = struct();
userdata.electric.egmX = [0.1 0.1 0; 0.8 0.1 0; 0.2 0.7 0];
userdata.electric.egmSurfX = userdata.electric.egmX;
userdata.electric.voltages = struct('bipolar', [1.1; 0.9; 1.4]);
end

function writeTextFile(filePath, text)
fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Could not create test file: %s', filePath);
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);
end

function tf = hasCheck(report, id)
tf = any(strcmp({report.checks.id}, id));
end
