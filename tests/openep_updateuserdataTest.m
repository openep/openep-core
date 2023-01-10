classdef openep_updateuserdataTest < matlab.unittest.TestCase

    properties(TestParameter)
       original_userdata = {'openep_dataset_1.mat'; load('openep_dataset_1.mat')};
       load_userdata = {true; false};
       version = {'-v6', '-v7', '-v7.3'};
       unpack_trirep = {{true, 'struct'}; {false, 'TriRep'}};
    end

    properties
        inputDir;
        outputPath;
    end

    methods(TestClassSetup)

        function determineInputDir(testCase)
            testCase.inputDir = fileparts(which('openep_dataset_1.mat'));
        end

        function determinOutputPath(testCase)
            import matlab.unittest.fixtures.TemporaryFolderFixture
            outputDir = testCase.applyFixture(TemporaryFolderFixture);
            testCase.outputPath = fullfile(outputDir.Folder, 'openep_dataset_1.mat');
        end

    end

    methods(Test, ParameterCombination = 'pairwise')

        function checkTriRepType(testCase, load_userdata, version, unpack_trirep)
            
            [unpack, triRepType] = unpack_trirep{:};
            
            if load_userdata
                load(fullfile(testCase.inputDir, 'openep_dataset_1.mat'), 'userdata');
            else
                userdata = fullfile(testCase.inputDir, 'openep_dataset_1.mat');
            end

            openep_updateuserdata( ...
                userdata, ...
                'version', version, ...
                'outputpath', testCase.outputPath, ...
                'unpacktrirep', unpack);
            
            new_userdata = load(testCase.outputPath, 'userdata');
            testCase.verifyTrue(isa(new_userdata.userdata.surface.triRep, triRepType));
        
        end

    end

end