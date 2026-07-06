# ARO/TRE Mapping Data Conversion

## Purpose

`convert_mapping_case` is the supported headless entry point for converting
CARTO and EnSiteX exports into OpenEP MAT files. It is intended for automated
execution where no MATLAB dialogs or progress windows are available.

The function:

1. prepares and validates the source export;
2. runs the appropriate importer without GUI interaction;
3. validates the imported OpenEP structure;
4. atomically publishes the MAT file;
5. writes JSON status and a text execution log.

An unsuccessful conversion does not replace an existing MAT output with a
partial file.

## CARTO

CARTO input may be an extracted export folder, its study XML, or a ZIP archive.
ZIP archives are inspected before extraction. On Linux, extraction uses
`/dev/shm` when sufficient storage and memory are available, then falls back to
the system temporary directory. Temporary files are removed after conversion.

Headless CARTO conversion requires explicit map, reference and ECG channels:

```matlab
addpath('/path/to/openep-core');

result = convert_mapping_case( ...
    '/input/carto-export.zip', ...
    '/output/case.mat', ...
    'system', 'carto', ...
    'maptoread', '2-LA', ...
    'refchannel', 'CS1-CS2', ...
    'ecgchannel', 'V1');
```

The MAT file contains a variable named `userdata`.

## EnSiteX

An EnSiteX study may contain separate exports for bipolar, unipolar and
omnipolar recordings. The converter discovers recording modes from file
contents and stores the selected modes as independent datasets:

```matlab
addpath('/path/to/openep-core');

result = convert_mapping_case( ...
    '/input/ensitex-study', ...
    '/output/case.mat', ...
    'system', 'ensitex', ...
    'maptoread', 'VoXel SR 1 ENDO', ...
    'modes', {'bi', 'uni', 'omni'});
```

The MAT file contains `openepCase`. Each entry in
`openepCase.datasets` contains one `userdata` structure and its recording mode,
source folder and detection evidence. Use:

```matlab
userdata = select_openep_dataset(openepCase, 'omni');
```

## Status Contract

By default, conversion produces three files:

```text
case.mat
case.status.json
case.log.txt
```

While conversion is running, it also maintains:

```text
case.progress.json
```

This file is atomically replaced as work advances. It contains the current
stage, percentage, message, update timestamp and elapsed time. CARTO reports
point, annotation, electrogram and force-data progress. EnSiteX reports
discovery and completion of each requested recording mode. The progress file
is removed only after the final status JSON has been written successfully. If
MATLAB is terminated unexpectedly, a stale progress file remains as evidence
of the interrupted job.

Progress can be watched from a shell:

```bash
watch -n 2 cat /output/case.progress.json
```

The JSON document is the machine-readable result. Important fields are:

| Field | Meaning |
|---|---|
| `success` | `true` only when import, output validation and MAT publication completed |
| `status` | `success`, `warning`, or `failure` |
| `outputFile` | Final MAT path |
| `outputPublished` | Whether this conversion published a new MAT file |
| `progressFile` | Path used for transient progress updates |
| `inputValidation` | Input checks and summary |
| `outputValidation` | OpenEP structure checks and summary |
| `runtimeWarning` | Last MATLAB importer warning, when present |
| `error` | Exception identifier, message and stack on failure |
| `timings` | Preparation, validation, import, save and total seconds |

All detailed validation checks have a stable `id`, a `level`, a stage, a
message and the related file. A successful conversion can have status
`warning`; the MAT output is still valid, but the warnings should be shown to
the user or retained for review.

The text log contains importer console output and warnings that are not part of
the structured validator. Paths and importer messages may contain study
identifiers, so status and log files must remain inside the TRE.

Custom status and log locations can be supplied with `statusfilename` and
`logfilename`.

## Batch Failure Behaviour

For a service that reads the JSON result, use the default
`throwonfailure=false`. The MATLAB call returns a result and writes diagnostics
even when conversion fails.

For a shell job that should exit nonzero:

```matlab
convert_mapping_case(inputPath, outputFile, ...
    'system', sourceSystem, ...
    'maptoread', mapName, ...
    'throwonfailure', true);
```

The status and log files are written before the exception is rethrown.

## Validation Levels

`validationlevel` accepts:

- `quick`: required files and lightweight format checks;
- `standard`: normal conversion checks and representative numeric validation;
- `full`: all currently implemented import-relevant checks.

Every converted OpenEP output receives structural and numeric validation
regardless of the input validation level.

## Verification

Run the normal unit and regression suite:

```matlab
run('run_validation_tests.m')
```

Run the complete CARTO integration test:

```bash
RUN_FULL_CARTO_IMPORT_TESTS=1 matlab -batch \
  "run('run_validation_tests.m')"
```

Run the complete EnSiteX multi-dataset integration test:

```bash
RUN_FULL_IMPORTER_SMOKE_TESTS=1 matlab -batch \
  "run('run_validation_tests.m')"
```

Full-case tests are opt-in because they require local clinical exports,
substantial storage and several minutes to complete.
