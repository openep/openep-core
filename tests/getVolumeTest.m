classdef getVolumeTest < matlab.unittest.TestCase
    methods(TestClassSetup)
        function loadData()
            load 'openep_dataset_2.mat' dataset_2;
        end
    end
    methods(Test)
        function dataset2Volume(testCase)
            expected_volume = 173.8810;
            calculated_volume = getVolume(dataset_2);
            testCase.verifyEqual(expected_volume,calculated_volume)
        end
    end
end