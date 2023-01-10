function newUserdata = doubleToSingle(userdata)
% DOUBLETOSINGLE stores data arrays as single rather than double
%
% Usage:
%   newUserdata = doubleToSingle(userdata)
% Where:
%   userdata - an OpenEP dataset, see https://openep.io/data/
%   newUserdata - a new OpenEP dataset
%
% DOUBLETOSINGLE converts the following arrays from double to single type:
%       userdata.electric.egmX
%       userdata.electric.egm
%       userdata.electric.egmUniX
%       userdata.electric.egmUni
%       userdata.electric.egmRef
%       userdata.electric.ecg
%       userdata.electric.annotations.woi
%       userdata.electric.annotations.referenceAnnot
%       userdata.electric.annotations.mapAnnot
%       userdata.electric.voltages.bipolar
%       userdata.electric.voltages.unipolar
%       userdata.electric.egmSurfX
%       userdata.electric.barDirection
%       userdata.electric.include
%
% Author: Steven Williams (2022) (Copyright)
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

% Copy userdata
newUserdata = userdata;

fieldToConvert = { 'electric.egmX' ...
    'electric.egm' ....
    'electric.egmUniX' ....
    'electric.egmUni' ....
    'electric.egmRef' ....
    'electric.ecg' ....
    'electric.annotations.woi' ....
    'electric.annotations.referenceAnnot' ....'' ....
    'electric.annotations.mapAnnot' ....
    'electric.voltages.bipolar' ....
    'electric.voltages.unipolar' ....
    'electric.egmSurfX' ....
    'electric.barDirection' ....
    'electric.include' ....
    };

for i = 1:numel(fieldToConvert)
    f = strsplit(fieldToConvert{i}, '.');
    switch(numel(f))
        case 1
            if isfield(userdata, f{1})
                newUserdata.(f{1}) = single(userdata.(f{1}));
            end
        case 2
            if isfield(userdata.(f{1}), f{2})
                newUserdata.(f{1}).(f{2}) = single(userdata.(f{1}).(f{2}));
            end
        case 3
            if isfield(userdata.(f{1}).(f{2}), f{3})
                newUserdata.(f{1}).(f{2}).(f{3}) = single(userdata.(f{1}).(f{2}).(f{3}));
            end
        case 4
            if isfield(userdata.(f{1}).(f{2}).(f{3}), f{4})
                newUserdata.(f{1}).(f{2}).(f{3}).(f{4}) = single(userdata.(f{1}).(f{2}).(f{3}).(f{4}));
            end
        otherwise
            error('OPENEP/DOUBLETOSINGLE: Not configured for more than 4 subfields');
    end
end
