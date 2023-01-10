classdef getVerticesTest < matlab.unittest.TestCase

    properties
        dataset_2
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.dataset_2 = load('openep_dataset_2.mat').userdata;
        end
    end

    methods(Test)
        function getVertices(testCase)

            v_test = getVertices(testCase.dataset_2);
            v_actual = testCase.dataset_2.surface.triRep.X;

            testCase.assertEqual(v_test, v_actual);

        end
    end
end