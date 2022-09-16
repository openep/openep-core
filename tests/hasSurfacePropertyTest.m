classdef hasSurfacePropertyTest < matlab.unittest.TestCase

    properties
        dataset_2
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.dataset_2 = load('openep_dataset_2.mat').userdata;
        end
    end

    methods(Test)
        function hasSurfaceProperty(testCase)

            % first we have to make sure a property exists
            mesh = getMesh(testCase.dataset_2, 'type', 'struct');
            numTriangles = size(mesh.Triangulation, 1);
            surfaceData = rand(numTriangles,1);
            propertyName = 'surfaceOfOrigin';
            propertyDefinedOn = 'elements';
            userdata = setSurfaceProperty(testCase.dataset_2 ...
                , 'name', propertyName ...
                , 'map', surfaceData ...
                , 'definedOn', propertyDefinedOn ...
                );

            % do the test
            tf = hasSurfaceProperty(userdata, propertyName);
            testCase.assertTrue(tf);

        end

        function doesNotHaveSurfaceProperty(testCase)
            propertyName = 'surfaceOfOrigin';

            tf = hasSurfaceProperty(testCase.dataset_2, propertyName);

            testCase.assertFalse(tf);

        end
    end
end