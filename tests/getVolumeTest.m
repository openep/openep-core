classdef getVolumeTest < matlab.unittest.TestCase

    properties
        dataset_2
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.dataset_2 = load('openep_dataset_2.mat').userdata;
        end
    end
    methods(Test)
        function dataset2Volume(testCase)
            expected_volume = 100;
            actual_volume = getVolume(testCase.dataset_2);
            testCase.verifyEqual(round(actual_volume, 3), expected_volume)
        end
    end
end