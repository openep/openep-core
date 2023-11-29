function newText = translateDataExportElement(oldText)
% TRANSLATEDATAEXPORTELEMENT standardises the description of Ensite files.
% The naming of files has changed, this causes problems. All filenames are
% mapped where possible.
% Usage:
%   newText = translateDataExportElement(oldText)
%   allDataExportTypes = translateDataExportElement()
%
% translator =   {
%            {'standardNameA', 'precisionName1A' , 'precisionName2A' , etc} ;
%            {'standardNameB', 'precisionName1B' , 'precisionName2B' , etc} ;
%               etc};
% allTranslations is a list of all {'standardNameA', 'standardNameB' , etc}

    newText = '';
    
    translator = { ...
        {'epcath_bip_raw' , 'EP_Catheter_Bipolar_Raw' , 'EPcathBIObipol_RAW' , 'EP_Catheter_Bipolar_Waveforms_Raw'};...
        {'epcath_bip_filt' , 'EP_Catheter_Bipolar_Filtered' , 'EPcathBIObipol_FILTERED' , 'EP_Catheter_Bipolar_Waveforms_Filtered'};...
        {'epcath_uni_raw' , 'EP_Catheter_Unipolar_Raw' , 'EPcathBIO_RAW' , 'EP_Catheter_Unipolar_Waveforms_Raw'};...
        {'epcath_uni_filt' , 'EP_Catheter_Unipolar_Filtered' , 'EPcathBIO_FILTERED' , 'EP_Catheter_Unipolar_Waveforms_Filtered'};...
        {'epcath_uni_comp' , 'EPcathBIO_COMPUTED'};...
        {'ecg_raw' , 'ECG_RAW' , 'ECG_Waveforms_Raw'}; ...
        {'ecg_filt' , 'ECG_FILTERED' , 'ECG_Waveforms_Filtered'}; ...
        {'respiration' , 'Respiration'}; ...
        {'locations' , 'Electrode_Locations' , 'Locations'}; ...
        {'channels' , 'Channels'}; ...
        {'dxldata' , '?????'}; ...
        {'modelgroups' , 'modelgroups' , 'Model_Groups'}; ...
                };
    
    if nargin == 0
        newText = cell(numel(translator),1);
        for i = 1:numel(translator)
            newText{i} = translator{i}{1};
        end
        return
    end

    for iTr=1:numel(translator)
        if any(contains(oldText,translator{iTr},'IgnoreCase',true))
            newText = translator{iTr}{1};
            break
        end
    end

end