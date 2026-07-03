classdef CartoFullCaseImportTest < matlab.unittest.TestCase
    % Opt-in integration test for a complete CARTO archive.

    methods (TestClassSetup)
        function addProjectPaths(~)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(repoRoot);
            addpath(fullfile(repoRoot, 'validation'));
        end
    end

    methods (Test)
        function importsStudy1Map2LAIntoValidatedUserdata(testCase)
            testCase.assumeTrue(runFullCartoTests(), ...
                'Set RUN_FULL_CARTO_IMPORT_TESTS=1 to run the full CARTO import.');

            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            casePath = getenv('OPENEP_FULL_CARTO_CASE');
            if isempty(casePath)
                casePath = fullfile(fileparts(repoRoot), 'full_cases', 'Carto', ...
                    'ReCETT-AF', 'Study1-Williams-Edinburgh', ...
                    'Export_PAF-02_28_2023-14-19-48.zip');
            end
            testCase.assumeTrue(isfolder(casePath) || isfile(casePath), ...
                'Full CARTO Study1 test case was not found.');

            mapName = environmentDefault('OPENEP_FULL_CARTO_MAP', '2-LA');
            refChannel = environmentDefault( ...
                'OPENEP_FULL_CARTO_REFCHANNEL', 'CS1-CS2');
            ecgChannel = environmentDefault( ...
                'OPENEP_FULL_CARTO_ECGCHANNEL', 'V1');

            preparationStart = tic;
            [caseFolder, extractionCleanup, extractionInfo] = ...
                prepare_carto_case_for_test(casePath);
            preparationSeconds = toc(preparationStart);
            studyXml = findStudyXml(caseFolder);

            importStart = tic;
            userdata = importcarto_mem(studyXml, ...
                'maptoread', mapName, ...
                'refchannel', refChannel, ...
                'ecgchannel', ecgChannel, ...
                'verbose', false);
            importSeconds = toc(importStart);

            validationStart = tic;
            report = validate_mapping_input(userdata, 'openep_userdata');
            validationSeconds = toc(validationStart);
            printTimings(extractionInfo, preparationSeconds, ...
                importSeconds, validationSeconds);

            testCase.verifyEqual(report.numFail, 0, report.summary);
            testCase.verifyTrue(isstruct(userdata.surface.triRep));
            verifyMesh(testCase, userdata.surface.triRep);
            verifyElectricData(testCase, userdata.electric, mapName);

            extractionRoot = extractionInfo.extractionRoot;
            delete(extractionCleanup);
            if extractionInfo.wasArchive
                testCase.verifyFalse(isfolder(extractionRoot));
            end
        end
    end
end

function tf = runFullCartoTests()
tf = any(strcmpi(getenv('RUN_FULL_CARTO_IMPORT_TESTS'), ...
    {'1', 'true', 'yes'}));
end

function value = environmentDefault(name, defaultValue)
value = getenv(name);
if isempty(value)
    value = defaultValue;
end
end

function studyXml = findStudyXml(caseFolder)
xmlFiles = dir(fullfile(caseFolder, '*.xml'));
names = {xmlFiles.name};
isStudyFile = ~startsWith(names, '.') & ...
    ~contains(names, 'Point_Export') & ...
    ~contains(names, 'Points_Export');
xmlFiles = xmlFiles(isStudyFile);

if numel(xmlFiles) ~= 1
    error('CartoFullCaseImportTest:StudyXmlSelection', ...
        'Expected one CARTO study XML in %s, found %d.', ...
        caseFolder, numel(xmlFiles));
end
studyXml = fullfile(xmlFiles(1).folder, xmlFiles(1).name);
end

function verifyMesh(testCase, mesh)
testCase.verifyTrue(isfield(mesh, 'X'));
testCase.verifyTrue(isfield(mesh, 'Triangulation'));
testCase.verifyGreaterThan(size(mesh.X, 1), 0);
testCase.verifyEqual(size(mesh.X, 2), 3);
testCase.verifyGreaterThan(size(mesh.Triangulation, 1), 0);
testCase.verifyEqual(size(mesh.Triangulation, 2), 3);
testCase.verifyTrue(all(isfinite(mesh.X), 'all'));

faces = mesh.Triangulation;
testCase.verifyTrue(all(isfinite(faces), 'all'));
testCase.verifyTrue(all(faces == round(faces), 'all'));
testCase.verifyGreaterThanOrEqual(min(faces, [], 'all'), 1);
testCase.verifyLessThanOrEqual(max(faces, [], 'all'), size(mesh.X, 1));
end

function verifyElectricData(testCase, electric, mapName)
nPoints = size(electric.egmX, 1);
if strcmp(mapName, '2-LA')
    testCase.verifyEqual(nPoints, 711);
else
    testCase.verifyGreaterThan(nPoints, 0);
end

testCase.verifySize(electric.egmX, [nPoints, 3]);
testCase.verifyEqual(size(electric.egm, 1), nPoints);
testCase.verifyEqual(size(electric.egmUni, 1), nPoints);
testCase.verifyEqual(size(electric.egmUni, 3), 2);
testCase.verifySize(electric.egmUniX, [nPoints, 3, 2]);
testCase.verifyEqual(size(electric.egmRef, 1), nPoints);
testCase.verifyEqual(size(electric.ecg, 1), nPoints);
testCase.verifyEqual(numel(electric.names), nPoints);
testCase.verifyEqual(size(electric.annotations.mapAnnot, 1), nPoints);
testCase.verifyEqual(size(electric.voltages.bipolar, 1), nPoints);
testCase.verifyEqual(size(electric.voltages.unipolar, 1), nPoints);
end

function printTimings(info, preparationSeconds, importSeconds, validationSeconds)
fprintf(['CARTO full-case timing: preparation %.1f s ', ...
    '(extraction %.1f s), import %.1f s, validation %.1f s.\n'], ...
    preparationSeconds, info.extractionSeconds, importSeconds, validationSeconds);
if info.wasArchive
    fprintf('CARTO archive: %d files, %.2f GiB uncompressed, root %s.\n', ...
        info.archiveFileCount, info.archiveUncompressedBytes / 1024^3, ...
        info.extractionRoot);
end
end
