classdef CVTriangulationTest < matlab.unittest.TestCase

    properties
        userdata
        cv
        cv_values
        cv_triangles
    end

    methods(TestClassSetup)
        function loadData(testCase)
            
            testCase.userdata = load('openep_dataset_2.mat').userdata;
            
            testCase.cv = CVTriangulation();
            testCase.cv.plot = false;
            [testCase.cv_values, testCase.cv_triangles] = testCase.cv.run(testCase.userdata);
            
        end
    end

    methods(Test)
        function nValues(testCase)
            
            expected_n_values = 669;
            testCase.verifyEqual(size(testCase.cv_values, 1), expected_n_values)

        end

        function firstCV(testCase)

            expected_first_cv = 0.1924;
            testCase.verifyEqual(round(testCase.cv_values(1), 4), expected_first_cv)

        end

        function firstTriangle(testCase)

            expected_first_triangle = [-13.2943, -101.7200,  118.9570];
            testCase.verifyEqual(testCase.cv_triangles(1, :), expected_first_triangle)

        end
    end
end
