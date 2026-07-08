function [caseFolder, cleanupObj, info] = prepare_carto_case_for_test(cartoPath)
%PREPARE_CARTO_CASE_FOR_TEST Compatibility wrapper for validation tests.

[caseFolder, cleanupObj, info] = prepare_carto_case(cartoPath);
end
