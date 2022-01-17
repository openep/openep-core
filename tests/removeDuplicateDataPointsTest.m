classdef removeDuplicateDataPointsTest < matlab.unittest.TestCase

    properties(TestParameter)
        x = {  [1 0.5 10; 5 6 8; 9 10 11; 1 0.5 10; 1 0.5 10] ...
            , [1 0.5 10; 5 6 8; 9 10 11; 1 0.5 10; 1 0.5 10] ...
            };
        f_x = { [5; 3; 4; 10; 12] ...
            , [6; 4; 3; 4; 2] ...
            };
        expected_x1 = { [1 0.5 10; 5 6 8; 9 10 11] ...
            , [1 0.5 10; 5 6 8; 9 10 11] ...
            };
        expected_f_x1 = { [12; 3; 4] ...
            , [4; 4; 3] ...
            };
        shouldWarn = {false ...
            , true ...
            }

        
    end

    methods(Test, ParameterCombination = 'sequential')
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

    %     methods(Test)
    %         x = [1 0.5 10; 5 6 8; 9 10 11; 1 0.5 10; 1 0.5 10];
    %         f_x = [6; 4; 3; 4; 2];
    %         expected
    %
    %         function testRemoveDuplicateDataPointsWithWarning(testCase, x, f_x)
    %         end
end