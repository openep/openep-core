classdef getVolumeTest < matlab.unittest.TestCase
    methods(Test)
        function dataset2Volume(testCase)
            load 'openep_dataset_2.mat' userdata;
            expected_volume = 173.8810;
            calculated_volume = getVolume(userdata);
            testCase.verifyEqual(expected_volume,calculated_volume)
        end
    end
end