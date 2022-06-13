classdef createUserDataTest < matlab.unittest.TestCase

    properties
        exampleUserData
    end

    methods(TestClassSetup)
        function setUp(testCase)

            d.systemName = [];
            d.cartoFolder = [];
            d.velocityFolder = [];
            d.precisionFolder = [];
            d.rhythmiaFolder = [];
            d.electric = [];
            d.notes = [];
            d.surface = [];
            d.rf = [];
            d.rfindex = [];

            d.electric.tags = [];
            d.electric.names = [];
            d.electric.electrodeNames_bip = [];
            d.electric.egmX = [];
            d.electric.egm = [];
            d.electric.electrodeNames_uni = [];
            d.electric.egmUniX = [];
            d.electric.egmUni = [];
            d.electric.egmRef = [];
            d.electric.ecg = [];
            d.electric.sampleFrequency = [];
            d.electric.annotations = [];
            d.electric.voltages = [];
            d.electric.impedances = [];
            d.electric.egmSurfX = [];
            d.electric.barDirection = [];
            d.electric.egmRefNames = [];
            d.electric.ecgNames = [];

            d.electric.annotations.woi = [];
            d.electric.annotations.referenceAnnot = [];
            d.electric.annotations.mapAnnot = [];

            d.electric.impedances.time = [];
            d.electric.impedances.value = [];

            d.electric.voltages.bipolar = [];
            d.electric.voltages.unipolar = [];

            d.surface.triRep = [];
            d.surface.isVertexAtRim = [];
            d.surface.act_bip = [];
            d.surface.uni_imp_frc = [];
            d.surface.signalMaps = [];
            d.surface.normals = [];
            d.surface.surfaceProperties = [];

            testCase.exampleUserData = d;
        end
    end

    methods(Test)
        function createUserData(testCase)

            newUserData = openep_createuserdata();

            isequaln(testCase.exampleUserData, newUserData)

            testCase.assertTrue(isequaln(testCase.exampleUserData, newUserData));

        end
    end
end