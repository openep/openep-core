classdef InspectEnsiteXExportTest < matlab.unittest.TestCase
    properties
        TempRoot
    end

    methods (TestClassSetup)
        function createTemporaryRoot(testCase)
            repoRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            addpath(repoRoot);
            testCase.TempRoot = tempname;
            mkdir(testCase.TempRoot);
        end
    end

    methods (TestClassTeardown)
        function removeTemporaryRoot(testCase)
            if isfolder(testCase.TempRoot)
                rmdir(testCase.TempRoot, 's');
            end
        end
    end

    methods (Test)
        function renamedFilesUseSemanticHeaders(testCase)
            folder = testCase.newExportFolder('renamed');
            writeMapCsv(fullfile(folder, 'arbitrary-a.csv'), 'LAT_bi', commonColumns());
            writeWaveCsv(fullfile(folder, 'arbitrary-b.csv'), 'refs');
            writeWaveCsv(fullfile(folder, 'arbitrary-c.csv'), 'rov');

            manifest = inspectensitex_export(testCase.TempRoot);
            export = manifest.exports(1);

            testCase.verifyEqual(export.recordingMode, 'bi');
            testCase.verifyEqual(export.confidence, 'high');
            testCase.verifyEqual(sort({export.waveFiles.role}), {'refs', 'rov'});
            testCase.verifyTrue(all(strcmp({export.waveFiles.source}, 'header')));
        end

        function waveFilenameFallbackIsReported(testCase)
            folder = testCase.newExportFolder('fallback');
            writeMapCsv(fullfile(folder, 'map.csv'), 'LAT_bi', commonColumns());
            writeWaveCsv(fullfile(folder, 'Wave_refs.csv'), '');

            manifest = inspectensitex_export(folder);

            testCase.verifyEqual(manifest.exports.waveFiles.role, 'refs');
            testCase.verifyEqual(manifest.exports.waveFiles.source, 'filename');
        end

        function omnipolarSchemaWorksWithoutModeSuffix(testCase)
            folder = testCase.newExportFolder('omni-schema');
            columns = [commonColumns(), {
                'pp_Vmax', 'pp_Valong', 'pp_Vacross', ...
                'Uni_Corner_Elec', 'Uni_Along_Elec', 'Uni_Across_Elec'
            }];
            writeMapCsv(fullfile(folder, 'map.csv'), 'LAT', columns);

            manifest = inspectensitex_export(folder);

            testCase.verifyEqual(manifest.exports.recordingMode, 'omni');
            testCase.verifyEqual(manifest.exports.confidence, 'high');
        end

        function insufficientEvidenceRemainsUnknown(testCase)
            folder = testCase.newExportFolder('unknown');
            writeMapCsv(fullfile(folder, 'map.csv'), 'LAT', commonColumns());
            writeWaveCsv(fullfile(folder, 'signal.csv'), 'rov');

            manifest = inspectensitex_export(folder);

            testCase.verifyEqual(manifest.exports.recordingMode, 'unknown');
        end

        function conflictingMapHeadersAreReported(testCase)
            folder = testCase.newExportFolder('conflict');
            writeMapCsv(fullfile(folder, 'map-one.csv'), 'LAT_bi', commonColumns());
            writeMapCsv(fullfile(folder, 'map-two.csv'), 'PP_uni', commonColumns());

            manifest = inspectensitex_export(folder);

            testCase.verifyEqual(manifest.exports.recordingMode, 'conflict');
            testCase.verifyNotEmpty(manifest.exports.errors);
        end
    end

    methods (Access = private)
        function folder = newExportFolder(testCase, name)
            folder = fullfile(testCase.TempRoot, name, 'anything');
            mkdir(folder);
            writeTextFile(fullfile(fileparts(folder), ...
                'Contact_Mapping_Model.xml'), '<model />');
        end
    end
end

function columns = commonColumns()
columns = {'Rov trace', 'Electrodes', 'Freeze Grp #', '(Point #)', ...
    'roving x', 'roving y', 'roving z', 'LAT'};
end

function writeMapCsv(filePath, mapType, columns)
header = {
    'Export File Version: 11'
    'Export Data Element: DxL'
    'Map name:,Test Map'
    ['Map type:,', mapType]
    '# mapping pts:,2'
};
writeDxlCsv(filePath, header, columns);
end

function writeWaveCsv(filePath, waveName)
header = {
    'Export File Version: 11'
    'Export Data Element: DxL'
    'Map name:,Test Map'
    ['Wave name:,', waveName]
    '# freeze groups:,2'
    'Sample rate:,2000'
};
writeDxlCsv(filePath, header, ...
    {'Trace', 'Freeze Grp #', '(Point #)', 'startTime (abs)', ...
    'rovTime (wave samples)', '0'});
end

function writeDxlCsv(filePath, header, columns)
dataStartRow = numel(header) + 3;
fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Could not create test CSV: %s', filePath);
cleanupObj = onCleanup(@() fclose(fid));
for i = 1:numel(header)
    fprintf(fid, '%s\n', header{i});
end
fprintf(fid, 'Data starts in row,%d\n', dataStartRow);
fprintf(fid, '*****\n');
fprintf(fid, '%s\n', strjoin(columns, ','));
delete(cleanupObj);
end

function writeTextFile(filePath, text)
fid = fopen(filePath, 'w');
assert(fid ~= -1, 'Could not create test file: %s', filePath);
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', text);
delete(cleanupObj);
end
