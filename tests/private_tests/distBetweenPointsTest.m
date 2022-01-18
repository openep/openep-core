classdef distBetweenPointsTest < matlab.unittest.TestCase

    properties (TestParameter)
        pointA = {[0, 0, 0], [0, 0, 0]};
        pointB = {...
            [5, 5, 5], ...
            [[5, 5, 5]; [10, 10, 10]], ...
            };
        expected_distance = {sqrt(75), [sqrt(75); sqrt(300)]};
        warn = {false, true};
    end
    
    methods(TestClassSetup)
        function changeDirToPrivate(testCase)
            cd('../../private/');
        end
    end

    methods(TestClassTeardown)
        function changeDirToTests(testCase)
            cd('../tests/private_tests/');
        end
    end

    methods(Test, ParameterCombination = 'sequential')

        function distances(testCase, pointA, pointB, expected_distance)
            
            actual_distance = distBetweenPoints(...
                pointA, ...
                pointB);
            testCase.verifyEqual( ...
                round(actual_distance, 3), ...
                round(expected_distance, 3));
        end

        function distancesWarn(testCase, pointA, pointB, warn)
            if warn
                verifyWarning(...
                    testCase, ...
                    @() distBetweenPoints(pointA, pointB, warn), ...
                    'OPENEP:distances');
            else
                verifyWarningFree(...
                    testCase, ...
                    @() distBetweenPoints(pointA, pointB, warn));
            end
        end
    
    end

    methods(Test)

        function unequalDimensionsError(testCase)
            verifyError( ...
                testCase, ...
                @() distBetweenPoints([0, 0, 0], [10, 10]), ...
                'OPENEP:distances');
        end

        function unequalNumPointsError(testCase)
            verifyError( ...
                testCase, ...
                @() distBetweenPoints([[0, 0, 0]; [0, 0, 0]], [10, 10]), ...
                'OPENEP:distances');
        end

    end

end