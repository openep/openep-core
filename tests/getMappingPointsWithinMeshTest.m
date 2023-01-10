classdef getMappingPointsWithinMeshTest < matlab.unittest.TestCase

    properties(TestParameter)
        tol = {0.01, 1000};
        expected_total_within = {1234, 1318};
    end

    properties
        userdata;
    end

    methods(TestClassSetup)
        function loadData(testCase)
            testCase.userdata = load('openep_dataset_1.mat').userdata;
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
        function tolerances(testCase, tol, expected_total_within)
            
            iPoint = getMappingPointsWithinMesh(...
                testCase.userdata, ...
                'tol', tol);

            actual_total_within = sum(iPoint);

            testCase.verifyEqual(actual_total_within, expected_total_within);
        end
    end
    
    methods(Test)
        function plotting(testCase, tol)

            % Turn off specific warning we know about
            
            % Warning: Some input points are not referenced by the triangulation
            warning('off', "MATLAB:triangulation:PtsNotInTriWarnId");

            % Warning: Error updating Light.
            % Exceedded the maximum number (8) of light sources.
            warning('off', "MATLAB:handle_graphics:exceptions:SceneNode");

            % Check the plot is produced without further warnings
            verifyWarningFree( ...
                testCase, ...
                @() getMappingPointsWithinMesh( ...
                    testCase.userdata, ...
                    'tol', tol, ...
                    'plot', true));
            
        end
    end

end