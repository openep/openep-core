function data = fixVoltageAnnotations(userdata)
% FIXVOLTAGEANNOTATIONS Fixes the uni/bip reversal identified in
% userdata.electric.voltages in August 2014. This function is not required
% for any userdata .mat files created after 27th August 2014.
%
% Usage:
%   userdata = fixVoltageAnnotations(userdata)
% Where:
%   userdata - is the output
%   userdata - is the input, or 'openfile'
%
% FIXVOLTAGEANNOTATIONS detailed description goes here.
%
% Author: Steven Williams (2014) (Copyright)
% SPDX-License-Identifier: Apache-2.0
%
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

filename = [];
pathname = [];
data = [];
if ischar(userdata)
    if strcmpi(userdata, 'openfile')
        [filename, pathname] = uigetfile('*.mat' ...
            , 'Pick one or more userdata .mat file(s)' ...
            , 'multiselect', 'on' ...
            );
        if isequal(filename,0) || isequal(pathname,0)
            disp('User pressed cancel')
            return
        else
            for i = 1:numel(filename)
                disp(['Loading ... ', fullfile(pathname, filename{i})])
                fullfilepath{i} = [pathname filename{i}];
                data{i} = load(fullfilepath{i}, 'userdata');
            end
        end
    else
        error('FIXVOLTAGEANNOTATIONS: Invalid input');
    end
else
    data{1} = userdata;
end

% Edited the files
for i = 1:numel(data)
    bipolar = data{i}.userdata.electric.voltages.unipolar;
    unipolar = data{i}.userdata.electric.voltages.bipolar;
    
    data{i}.userdata.electric.voltages.bipolar = bipolar;
    data{i}.userdata.electric.voltages.unipolar = unipolar;
    if isfield(data{i}.userdata, 'notes')
        data{i}.userdata.notes{end+1} = [date ': ran fixVoltageAnnotations.m'];
    else
        data{i}.userdata.notes{1} = [date ': ran fixVoltageAnnotations.m'];
    end
end

% Save the files
for i = 1:numel(data)
    userdata = data{i}.userdata;
    save(fullfilepath{i}, 'userdata');
    disp(['Saved ...', fullfilepath{i}]);
end