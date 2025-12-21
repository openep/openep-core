classdef loadensitex_dxldata_SmokeTest < matlab.unittest.TestCase
    % loadensitex_dxldata_SmokeTest
    %
    % Smoke test for loadensitex_dxldata.
    %
    % This test verifies that all EnsiteX DXL CSV files in the OpenEP test
    % dataset can be parsed without throwing errors. The intent is not to
    % validate data correctness, but to ensure that the importer is robust
    % to a variety of real-world input files.
    %
    % Key behaviour:
    %   - Iterates over all subfolders in EnsiteXDXLFiles
    %   - Attempts to load every non-hidden *.csv file
    %   - Ignores warnings (errors only cause failures)
    %   - Treats empty dataset folders as failures
    %   - Optionally runs in parallel (requires Parallel Computing Toolbox)
    %   - Disables progress UI for headless/parallel execution
    %   - Reports a full pass/fail summary with timing information
    %

    properties (Constant)
        useParallel = true;   % Set true to enable parallel execution
    end

    methods (Test)
        function runsOnAllDXLDatasets(testCase)

            % Start overall timer
            tTotal = tic;

            % Locate test data
            testDataRoot = getenv("OPENEP_TESTING_DATA");
            testCase.assumeNotEmpty(testDataRoot, ...
                "OPENEP_TESTING_DATA environment variable not set.");

            rootDir = fullfile(testDataRoot, 'EnsiteXDXLFiles');
            testCase.assumeTrue(isfolder(rootDir), ...
                "EnsiteXDXLFiles directory not found.");

            % Find dataset folders
            datasets = dir(rootDir);
            datasets = datasets([datasets.isdir]);
            datasets = datasets(~ismember({datasets.name}, {'.','..'}));
            testCase.assumeNotEmpty(datasets, ...
                "No EnsiteX DXL datasets found.");

            % Collect CSV files and empty folders
            allCSVFiles  = {};
            emptyFolders = {};

            for i = 1:numel(datasets)
                dxlDir = fullfile(rootDir, datasets(i).name);
                csvFiles = dir(fullfile(dxlDir, '*.csv'));

                % Ignore hidden / macOS metadata files
                csvFiles = csvFiles(~startsWith({csvFiles.name}, '.') & ...
                                    ~startsWith({csvFiles.name}, '._'));

                if isempty(csvFiles)
                    emptyFolders{end+1} = datasets(i).name; %#ok<AGROW>
                else
                    for j = 1:numel(csvFiles)
                        allCSVFiles{end+1} = fullfile(datasets(i).name, csvFiles(j).name); %#ok<AGROW>
                    end
                end
            end

            numFiles = numel(allCSVFiles);
            passedFiles = cell(1, numFiles);
            failedFiles = cell(1, numFiles);
            durations   = zeros(1, numFiles);

            % -----------------------------------------------------------------
            % Load CSV files (parallel or sequential)
            % -----------------------------------------------------------------
            if testCase.useParallel
                parfor k = 1:numFiles
                    csvPath = fullfile(rootDir, allCSVFiles{k});

                    tStart_local = tic;

                    try
                        loadensitex_dxldata(csvPath, 'ShowProgress', false);
                        durations(k) = toc(tStart_local);
                        passedFiles{k} = struct( ...
                            'file', allCSVFiles{k}, ...
                            'duration', durations(k));
                    catch ME
                        durations(k) = toc(tStart_local);
                        failedFiles{k} = struct( ...
                            'file', allCSVFiles{k}, ...
                            'error', ME.message, ...
                            'duration', durations(k));
                    end
                end
            else
                for k = 1:numFiles
                    csvPath = fullfile(rootDir, allCSVFiles{k});
                    try
                        tStart_local = tic;
                        loadensitex_dxldata(csvPath, 'ShowProgress', false);
                        durations(k) = toc(tStart_local);
                        passedFiles{k} = struct( ...
                            'file', allCSVFiles{k}, ...
                            'duration', durations(k));
                    catch ME
                        durations(k) = toc(tStart_local);
                        failedFiles{k} = struct( ...
                            'file', allCSVFiles{k}, ...
                            'error', ME.message, ...
                            'duration', durations(k));
                    end
                end
            end

            % Clean up results
            passedFiles = passedFiles(~cellfun('isempty', passedFiles));
            failedFiles = failedFiles(~cellfun('isempty', failedFiles));

            % Add empty-folder failures
            for k = 1:numel(emptyFolders)
                failedFiles{end+1} = struct( ...
                    'file', emptyFolders{k}, ...
                    'error', 'No CSV files found in folder (after filtering hidden files)', ...
                    'duration', NaN); %#ok<AGROW>
            end

            % -----------------------------------------------------------------
            % Summary output
            % -----------------------------------------------------------------
            fprintf('\n================================================================================\n')
            fprintf('\n=== loadensitex_dxldata Smoke Test Summary ===\n');

            if ~isempty(passedFiles)
                fprintf('  Passed files (%d):\n', numel(passedFiles));
                for k = 1:numel(passedFiles)
                    fprintf('    %s  | Time: %.2f s\n', ...
                        passedFiles{k}.file, passedFiles{k}.duration);
                end
            else
                fprintf('  No CSV files loaded successfully.\n');
            end

            if ~isempty(failedFiles)
                fprintf('  Failed files (%d):\n', numel(failedFiles));
                for k = 1:numel(failedFiles)
                    if isnan(failedFiles{k}.duration)
                        durStr = 'N/A';
                    else
                        durStr = sprintf('%.2f s', failedFiles{k}.duration);
                    end
                    fprintf('    %s  | Error: %s | Time: %s\n', ...
                        failedFiles{k}.file, failedFiles{k}.error, durStr);
                end

                totalTime = toc(tTotal);
                fprintf('Total test runtime: %.2f seconds\n', totalTime);
                fprintf('=== End of Summary ===\n\n');
                fprintf('\n================================================================================\n')

                testCase.assertFail( ...
                    "One or more CSV files or folders failed. See summary above.");
            else
                fprintf('=== All CSV files loaded successfully ===\n');
                totalTime = toc(tTotal);
                fprintf('Total test runtime: %.2f seconds\n', totalTime);
                fprintf('=== End of Summary ===\n\n');
                fprintf('\n================================================================================\n')
            end

        end
    end
end
