function userdata = openep_createuserdata()
% OPENEP_CREATEUSERDATA creates the userdata for OpenEP.
%
% Usage:
%   userdata = openep_createuserdata(userdata)
% Output:
%   info - see 'format' below for details
% Inputs:
%
% OPENEP_CREATEUSERDATA creates an empty userdata structure
%
% Author: Steven Williams, Nick Linton (2020)
%
% Info on Code Testing:
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

    fNames = {    'systemName','cartoFolder','velocityFolder','precisionFolder' ...
                ,'rhythmiaFolder','electric','notes','surface','rf','rfindex'};
    for i = 1:numel(fNames)
        userdata.(fNames{i}) = [];
    end

    fNames = {   'tags','names','electrodeNames_bip','egmX','egm','electrodeNames_uni'...
                ,'egmUniX','egmUni','egmRef','ecg','sampleFrequency','annotations','voltages'...
                ,'impedances','egmSurfX','barDirection'};

    for i = 1:numel(fNames)
        tempstruct.(fNames{i}) = [];
    end

    annotationsSubFields = {'woi','referenceAnnot','mapAnnot'};
    voltagesSubFields = {'bipolar','unipolar'};
    impedancesSubFields = {'time','value'};

    for i = 1:numel(annotationsSubFields)
        tempstruct.annotations.(annotationsSubFields{i}) = [];
    end
    for i = 1:numel(voltagesSubFields)
        tempstruct.voltages.(voltagesSubFields{i}) = [];
    end
    for i = 1:numel(impedancesSubFields)
        tempstruct.impedances.(impedancesSubFields{i}) = [];
    end


    userdata.electric = tempstruct;


    fNames = {   'triRep','isVertexAtRim','act_bip','uni_imp_frc','signalMaps'};
    for i = 1:numel(fNames)
        tempstruct.(fNames{i}) = [];
    end
    userdata.surface = tempstruct;

end
