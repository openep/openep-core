classdef importPrecisionTest < matlab.unittest.TestCase

    properties
        precision_study;
    end
    
    properties (TestParameter)
        fields = {'precisionFolder', 'surface', 'electric', 'rf', 'notes'};
        surface_fields = {"triRep", "isVertexAtRim", "act_bip", "uni_imp_frc"};
        electric_fields = {...
            "tags", "names", "electrodeNames_bip", "egmX", "egm", ...
            "electrodeNames_uni", "egmUniX", "egmUni", "egmRef", "ecg", ...
            "annotations", "voltages", "impedances", ...
            "egmSurfX", "barDirection" ...
            };
        % rf_fields = {"force", "ablparams"};
    end

    methods (TestClassSetup)

        function loadPrecisionStudy(testCase)
            testCase.precision_study = importPrecisionTest.loadPrecision( ...
                getenv('PRECISION_TEST_CASE_1') ...
                );
        end

    end

    methods (Static)

        function precision = loadPrecision(study_file)
            study_file = fullfile(getenv("OPENEP_TESTING_DATA"), study_file);
            precision = import_precision( ... % this is the function we are testing
                study_file ...
                );
        end

    end

    methods (Test)

        function hasFields(testCase, fields)
            testCase.assertTrue(isfield(testCase.precision_study, fields));
        end

        function hasSurfaceFields(testCase, surface_fields)
            testCase.assertTrue(isfield(testCase.precision_study.surface, surface_fields));
        end

        function hasElectricFields(testCase, electric_fields)
            testCase.assertTrue(isfield(testCase.precision_study.electric, electric_fields));
        end

%         function hasRFFields(testCase, rf_fields)
%             testCase.assertTrue(isfield(testCase.precision_study.rf.originaldata, rf_fields))
%         end

        function dimensionTriRep(testCase)

            n_points = 106586;
            n_triangles = 213172;

            testCase.verifyEqual( ...
                [n_points, 3], ...
                size(testCase.precision_study.surface.triRep.X) ...
                );
            testCase.verifyEqual( ...
                [n_triangles, 3], ...
                size(testCase.precision_study.surface.triRep.Triangulation) ...
                );

        end

    end

end