classdef batchConvertTest < matlab.unittest.TestCase

    properties
        inputDir;
        outputDir;
        allOriginalFiles;
        allNewFiles;
    end

    methods(TestClassSetup)

        function determineInputDir(testCase)
            testCase.inputDir = fileparts(which('openep_dataset_1.mat'));
        end

        function determinOutputDir(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            testCase.outputDir = testCase.applyFixture(TemporaryFolderFixture);
        end

        function runBatchConvert(testCase)
            batchConvert(testCase.inputDir, testCase.outputDir.Folder);
        end

        function getFilenames(testCase)

            % We can't access openep-core/private/nameFile.m from this
            % folder
            cwd = pwd;

            filePath = which(mfilename);
            testFolder = fileparts(filePath);
            privateFolder = fullfile(testFolder, '../private');
            cd(privateFolder);

            testCase.allOriginalFiles = nameFiles( ...
                testCase.inputDir, ...
                'showhiddenfiles', false, ...
                'extension', 'mat');

            testCase.allNewFiles = nameFiles( ...
                testCase.outputDir.Folder, ...
                'showhiddenfiles', false, ...
                'extension', 'mat');

            cd(cwd);
        end

    end

    methods(Test)

        function checkFileNames(testCase)
            
            for i = 1:numel(testCase.allNewFiles)
                testCase.verifyEqual(testCase.allNewFiles{i}, testCase.allOriginalFiles{i});
            end
        
        end

        function checkTriRepType(testCase)
            
            for i = 1:numel(testCase.allNewFiles)
                filename = fullfile(testCase.outputDir.Folder, testCase.allNewFiles{i});
                load(filename, 'userdata');
                testCase.verifyTrue(isstruct(userdata.surface.triRep));
            end
        
        end

    end

end
