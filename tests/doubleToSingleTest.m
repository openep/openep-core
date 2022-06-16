classdef doubleToSingleTest < matlab.unittest.TestCase

    properties
        dataset_2
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.dataset_2 = load('openep_dataset_2.mat').userdata;
        end
    end

    methods(Test)
        function double2Single(testCase)

            userdata = doubleToSingle(testCase.dataset_2);

            testCase.assertTrue(isa(userdata.electric.egmX, "single"));
            
        end
    end
end