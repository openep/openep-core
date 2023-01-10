classdef getSurfacePropertyTest < matlab.unittest.TestCase

    properties
        dataset_2
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.dataset_2 = load('openep_dataset_2.mat').userdata;
        end
    end

    methods(Test)
        function getSurfaceProperty(testCase)

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
            calcSurfaceData = getSurfaceProperty(userdata, propertyName);
            testCase.assertEqual(calcSurfaceData.name, propertyName);
            testCase.assertEqual(calcSurfaceData.map, surfaceData);
            testCase.assertEqual(calcSurfaceData.definedOn, propertyDefinedOn);

        end
    end
end