classdef importCartoMemTest < matlab.unittest.TestCase

    properties
        carto_study1_map2;
    end
    
    properties (TestParameter)
        fields = {'cartoFolder', 'surface', 'electric', 'rf', 'notes'};
        surface_fields = {"triRep", "isVertexAtRim", "act_bip", "uni_imp_frc"};
        electric_fields = {...
            "tags", "names", "electrodeNames_bip", "egmX", "egm", ...
            "electrodeNames_uni", "egmUniX", "egmUni", "egmRef", "ecg", ...
            "annotations", "voltages", "impedances", ...
            "egmSurfX", "barDirection" ...
            };
        rf_fields = {"force", "ablparams"};
    end

    methods (TestClassSetup)

        function loadCartoStudy1Map2(testCase)
            testCase.carto_study1_map2 = importCartoMemTest.loadCarto( ...
                'Carto/Export_Study-1-11_25_2021-15-01-32/Study 1 11_25_2021 15-01-32.xml', ...
                '2-Map', ...
                'v1', ...
                'v1' ...
                );
        end

    end

    methods (Static)

        function carto = loadCarto(study_file, map, refchannel, ecgchannel)
            study_file = fullfile(getenv("OPENEP_TESTING_DATA"), study_file)
            carto = importcarto_mem( ...
                study_file, ...
                'maptoread', map, ...
                'refchannel', refchannel, ...
                'ecgchannel', ecgchannel, ...
                'savefilename', fullfile(getenv("OPENEP_TESTING_DATA"), 'tmp') ...
                );
        end

    end

    methods (Test)

        function hasFields(testCase, fields)
            testCase.assertTrue(isfield(testCase.carto_study1_map2, fields));
        end

        function hasSurfaceFields(testCase, surface_fields)
            testCase.assertTrue(isfield(testCase.carto_study1_map2.surface, surface_fields));
        end

        function hasElectricFields(testCase, electric_fields)
            testCase.assertTrue(isfield(testCase.carto_study1_map2.electric, electric_fields));
        end

        function hasRFFields(testCase, rf_fields)
            testCase.assertTrue(isfield(testCase.carto_study1_map2.rf.originaldata, rf_fields))
        end

        function dimensionTriRep(testCase)

            n_points = 2916;
            n_triangles = 5828;

            testCase.verifyEqual( ...
                [n_points, 3], ...
                size(testCase.carto_study1_map2.surface.triRep.X) ...
                );
            testCase.verifyEqual( ...
                [n_triangles, 3], ...
                size(testCase.carto_study1_map2.surface.triRep.Triangulation) ...
                );

        end

    end

end
