classdef removeDuplicateDataPointsTest < matlab.unittest.TestCase

    properties(TestParameter)
        % Set up input and expected output data
        x = {  [1 0.5 10; 5 6 8; 9 10 11; 1 0.5 10; 1 0.5 10] ... % test a duplicate which has a maximum value
            , [1 0.5 10; 5 6 8; 9 10 11; 1 0.5 10; 1 0.5 10] ... % test a duplicate which has both a maximum and minium value; and therefore gives a warning
            , [1 0.5 10; 5 6 8; 9 10 11; 1 0.5 10; 1 0.5 10] ... % test a minimum value
            , [1 0.5 10; 5 6 8; 9 10 11; 1 0.5 10; 1 0.5 10] ... % test an average
            };
        f_x = { [5; 3; 4; 10; 12] ... % test a duplicate which has a maximum value
            , [6; 4; 3; 4; 2] ... % test a duplicate which has both a maximum and minium value; and therefore gives a warning
            , [5; 8; 4; 2; 1] ... % test a minimum value
            , [6; 1; 10; 4; 2] ... % test an average
            };
        expected_x1 = { [1 0.5 10; 5 6 8; 9 10 11] ...
            , [1 0.5 10; 5 6 8; 9 10 11] ...
            , [1 0.5 10; 5 6 8; 9 10 11] ...
            , [1 0.5 10; 5 6 8; 9 10 11] ...
            };
        expected_f_x1 = { [12; 3; 4] ...
            , [4; 4; 3] ...
            , [1; 8; 4] ...
            , [4; 1; 10] ...
            };
        shouldWarn = {false ...
            , true ...
            , false ...
            , false ...
            }
    end

    methods(Test, ParameterCombination = 'sequential')
        % perform tests
        function testRemoveDuplicateDataPoints(testCase, x, f_x, expected_x1, expected_f_x1, shouldWarn)
            [x1, f_x1] = removeDuplicateDataPoints(x, f_x);
            testCase.verifyEqual(x1, expected_x1);
            testCase.verifyEqual(f_x1, expected_f_x1);

            if ~shouldWarn
                verifyWarningFree(testCase, @() removeDuplicateDataPoints(x,f_x));
            else
                verifyWarning(testCase, @() removeDuplicateDataPoints(x,f_x), 'OPENEP:information');
            end
        end
    end
end