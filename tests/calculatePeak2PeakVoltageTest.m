classdef calculatePeak2PeakVoltageTest < matlab.unittest.TestCase
    methods (Test)
        function calculatesVoltageWithoutOptionalToolboxes(testCase)
            egms = [
                0 1 5 2 0
                0 4 4 4 0
            ];
            referenceAnnotations = [3; 3];
            windowsOfInterest = [-1 1; -1 1];

            voltage = calculatePeak2PeakVoltage( ...
                egms, referenceAnnotations, windowsOfInterest);

            testCase.verifyEqual(voltage, [4; 0]);
        end
    end
end
