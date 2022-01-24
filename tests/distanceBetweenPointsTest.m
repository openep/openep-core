classdef distanceBetweenPointsTest < matlab.unittest.TestCase

    properties
        dataset_1
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.dataset_1 = load('openep_dataset_1.mat').userdata;
        end
    end
    methods(Test)
        function calculateLinearDistance(testCase)
            expected_distance = 28.527;
            actual_distance = distanceBetweenPoints(testCase.dataset_1, 1, 200, 'method', 'linear', 'plot', false);
            testCase.verifyEqual(round(actual_distance, 3), expected_distance)
        end
        function calculateGeodesicDistance(testCase)
            expected_distance = 32.610;
            actual_distance = distanceBetweenPoints(testCase.dataset_1, 1, 200, 'method', 'geodesic', 'plot', false);
            testCase.verifyEqual(round(actual_distance, 3), expected_distance)
        end
    end
end