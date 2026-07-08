function [userdata, dataset] = select_openep_dataset(openepCase, selector)
%SELECT_OPENEP_DATASET Select one userdata dataset by ID or recording mode.

assert(isstruct(openepCase) && isfield(openepCase, 'datasets') && ...
    isstruct(openepCase.datasets), ...
    'SELECT_OPENEP_DATASET: Invalid OpenEP case container.');

selector = char(selector);
ids = {openepCase.datasets.id};
modes = {openepCase.datasets.recordingMode};
matches = strcmpi(ids, selector) | strcmpi(modes, selector);
if sum(matches) ~= 1
    error('SELECT_OPENEP_DATASET: Selector must identify exactly one dataset.');
end

dataset = openepCase.datasets(matches);
userdata = dataset.userdata;
end
