classdef distanceBetweenPointsTest < matlab.unittest.TestCase

    properties(TestParameter)
        method = {'linear', 'geodesic'};
        expected_distance = {25.562, 27.4};
    end

    properties
        dataset_1;
        pointA;
        pointB;
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.dataset_1 = load('openep_dataset_1.mat').userdata;
            testCase.pointA = 1;
            testCase.pointB = 100;
        end
    end

    % Don't show plots when testing. Snippet from:
    % https://uk.mathworks.com/matlabcentral/answers/124174-suppress-figures-in-unit-tests#answer_131776
    methods(TestMethodSetup)
        function hidePlots(testCase)
            origDefault = get(0,'DefaultFigureVisible');
            set(0,'DefaultFigureVisible','off');
            testCase.addTeardown(@set, 0, 'DefaultFigureVisible', origDefault);
        end
    end

    methods(Test, ParameterCombination = 'sequential')
        function dataset1Distances(testCase, method, expected_distance)
            
            actual_distance = distanceBetweenPoints( ...
                testCase.dataset_1, ...
                testCase.pointA, ...
                testCase.pointB, ...
                'method', method);

            testCase.verifyEqual(round(actual_distance, 3), expected_distance);
        end
    end
    
    methods(Test)
        function dataset1Plot(testCase, method)
            verifyWarningFree( ...
                testCase, ...
                @() distanceBetweenPoints( ...
                    testCase.dataset_1, ...
                    testCase.pointA, ...
                    testCase.pointB, ...
                    'method', method, ...
                    'plot', true));
        end
    end

end