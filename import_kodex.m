function [userdata] = import_kodex(varargin)
% import_kodex provides a openEP data structure from multiple kodex files
% Usage:
%   import_kodex
%   userdata=import_kodex()
%   userdata=import_kodex('study_dir_overall')
%   userdata=import_kodex('study_dir_overall', optional parameter pairs.....)
% Where
%  study_dir_overall is the direcotry string with all files corresponding to an exported kodex study ()
%  (study_dir_overall can include multiple chambers/studies that can be chosen by parameter-pair input or pop up user input)
%  if no input, i.e. import_kodex or userdata=import_kodex(), a user
%  interafce will pop up to chose overall study directory,
%  userdata is a single openEP data structure
%
% import_kodex imports kodex data of formats:
%   1) Folder -> dump state -> main_entry -> studies -> 'LAT' with corresponding 'report' folder
%   2) Folder -> studies -> 'LAT'
%   For case 1, 'LAT' contains json files, while 'report folder contains csv points file (lat_points)
%   and mesh file (.obj) ('full' data output from kodex)
%   For case 2, 'LAT' contains csv points file (lat_points), mesh file (.obj) and (optional) .json files with raw electrogram data
%
% import_kodex accepts the following parameter-value pairs
%   'maptoread'     {''}|string
%       Specifies which map to read.
%  'savefilename'       {''}|string
%       The full path to the location in which to save the output.
%  'save_option'  true|{false}
%       Option for save dialog to pop up at the end. If a value is given for 'savefilename',
%       two versions of the file will be saved
%
% userdata structure ... (needs to be updated)
%   .surface
%       .triRep         - TriRep object for the surface
%       .isVertexAtRim  - logical array indicating vertices at a 'rim'
%       .act_bip        - nVertices*2 array of activation and voltage data
%       .uni_imp_frc    - nVertices*3 array of uni voltage, impedance and contact force
%   .electric
%       .isPointLocationOnly    - logical array
%       .tags
%       .names
%       .egmX           - location of point
%       .egmSurfX       - location of surface nearest point
%       .barDirection   - normal to surface at egmSurfX
%       .egm            - bipolar electrogram
%       .egmUni         - matrix of unipolar electrograms
%       .egmUniX        - location of unipolar points
%       .egmRefNames    - names of egmRef
%       .egmRef         - electrogram of reference
%       .ecgNames       - ecg names (or other channel names)
%       .ecg            - ecg
%       .force
%           .force    - instantaneous force recording
%           .axialAngle    - axial angle
%           .lateralAngle  - lateral angle
%           .time_force - time course of force [(:,:,1)=time, (:,:,2)=force]
%           .time_axial - time course of axial angle [(:,:,1)=time, (:,:,2)=axial angle]
%           .time_lateral - time course of lateral angle [(:,:,1)=time, (:,:,2)=lateral angle]
% NB: Frequency not yet read in
%% Upload

if nargin == 0
    study_dir_overall=uigetdir;
else
    study_dir_overall=varargin{1};
end
chamber={};
savefilename={};
save_option=false;
nStandardArgs=1;
if nargin > nStandardArgs +1
    for i = 2:2:nargin
        switch varargin{i}
            case 'maptoread'
                chamber=varargin{i+1};
            case 'savefilename'
                save_path=varargin{i+1};
            case 'save_option'
                save_option=varargin{i+1};
        end
    end
end


tempfolder=dir(study_dir_overall);

studies=[];
map_count=0;
for j=3:numel(tempfolder)
    if tempfolder(j).isdir == 1
        map_count=map_count+1;
        studies{map_count}=tempfolder(j).name;
    end
end

if map_count==0
    warning('No studies found')
end

%check format and if study directory needs to be changed to 'full' format
%could be removed if only expecting 'full' data exports
if any(strcmp(studies,'dump_state')) && any(strcmp(studies,'Report'))
    study_dir_new=[study_dir_overall,filesep(),'dump_state',filesep(),'MainEntry',filesep()]; %update study_dir
    tempfolder=dir(study_dir_new);
    map_count=0;
    studies=[];
    for j=3:numel(tempfolder)
        if tempfolder(j).isdir == 1
            map_count=map_count+1;
            studies{map_count}=tempfolder(j).name;
        end
    end
end

if exist('study_dir_new') == 0
    if isempty(chamber) == 1
        [indx,tf] = listdlg('Liststring',studies);
        chamber=char(studies(indx));
    end
    egm_folder=[study_dir_overall,filesep(),chamber,filesep(),'LAT'];
    abl_folder=[study_dir_new,filesep(),chamber,filesep(),'ABL'];
    lm_folder=[study_dir_new,filesep(),chamber,filesep(),'LM'];
    tempfolder=dir(egm_folder);
    tempfolder_abl=dir(abl_folder);
    tempfolder_lm=dir(lm_folder);
    mesh_folder=[study_dir_overall,filesep(),chamber,filesep(),'LAT'];
    files=dir(mesh_folder);
elseif exist('study_dir_new') == 1
    if isempty(chamber) == 1
        [indx,tf] = listdlg('Liststring',studies);
        chamber=char(studies(indx))
    end
    egm_folder=[study_dir_new,filesep(),chamber,filesep(),'LAT'];
    abl_folder=[study_dir_new,filesep(),chamber,filesep(),'ABL'];
    lm_folder=[study_dir_new,filesep(),chamber,filesep(),'LM'];
    tempfolder=dir(egm_folder);
    tempfolder_abl=dir(abl_folder);
    tempfolder_lm=dir(lm_folder);
    mesh_folder=[study_dir_overall,filesep(),'Report',filesep(),chamber];
    files=dir(mesh_folder);
end


% Fs=500; %work out way to read frequency
% sampling_rate=1000/Fs; %in ms
userdata = openep_createuserdata();
%% Get software version and other details %%TODO: Make getter function?
userdata.kodexFolder=study_dir_overall;
task_folder=[study_dir_overall,filesep(),'task_manager_logs'];
task_dir=dir(task_folder);
task_folder=[task_folder,filesep(),task_dir(3).name];
task_file=[task_folder,filesep(),'MainNoPatientLoader.txt'];
fid = fopen(task_file);
tline = fgetl(fid);
lineCounter = 1;
while ischar(tline)
    if contains(tline,'system_label')
        break;
    end
    % Read next line
    tline = fgetl(fid);
    lineCounter = lineCounter + 1;
end
fclose(fid);
k = strfind(tline,'system_label');
kodex_version=tline(71:84);
userdata.kodexData.softwareVersion=kodex_version;

%% find object and report data files (as measured by kodex)
for n=3:numel(files)
    fname=files(n).name;
    [a,b,c]=fileparts(fname);
    if strcmp('.obj',c) == 1
        object_file=[mesh_folder,filesep(),b,c];
    end
    if strcmp('.csv',c) == 1 && strcmp('lat_points',b) == 1
        csv_file=[mesh_folder,filesep(),b,c];
    end
    if strcmp('.csv',c) == 1 && strcmp('landmark_points',b) == 1
        csv_landmark_file=[mesh_folder,filesep(),b,c];
    end
end



%% Get data and mesh from kodex analysis
T = readtable(csv_file);
W = width(T);
Values=T(:,W);
if size(Values,1) == 0
    warning('No data contained in Map');
end
datas=split(table2array(Values),',');
LATs=str2double(datas(:,5));
LATs=LATs-min(LATs);
Volts=str2double(datas(:,6));
csv_pos=str2double(datas(:,2:4));
obj = readObj(object_file);


% maybe osbsolte in 'full' files as data save in json files
T1=readtable(csv_landmark_file);
W1 = width(T1);
H1 = height(T1);
Landmarks=T1(:,W1);
if size(Landmarks,1) == 0
    warning('No landmark data contained in Map');
    Lpoints=[];
else
%datas2=split(table2array(Landmarks),',')
%Lpoints=str2double(datas2(:,2:4));
end
%% Get 'raw' electirc data
egm_count=0;
abl_count=0;
lm_count=0;

%% this is to fill in from 'raw' electrogram and ablation data if avalable (think of better way)


%userdata.electric.sampleFrequency=0; %TO DO - Find where this is stored

if exist('study_dir_new') == 0

%egm data    
for i=1:numel(tempfolder);
    fname=tempfolder(i).name;
    [a,b,c]=fileparts(fname);
    if strcmp(c,'.json') == 1 && contains(b,'LAT')
                egm_count=egm_count+1;
                x=kodex_egm_import([egm_folder,filesep(),fname]);
                value = jsondecode(x);
                
                index(egm_count)=value.index_on_mesh;
                
                userdata.electric.indexmesh(egm_count,1)=value.index_on_mesh;
                
                %electric data
                userdata.electric.tags{egm_count,1}=value.tag;
                userdata.electric.names{egm_count,1}=value.current_point;
                userdata.electric.egmX(egm_count,:)=value.position; %%is raw position
                userdata.electric.voltages.bipolar(egm_count,1)=Volts(value.current_point); %NB-Unipolar data not currently collected by Kodex
                userdata.electric.LATs(egm_count,1)=LATs(value.current_point); %keep for the moment but check correlation with getting LAT from annotations, if they match can remove
                userdata.electric.egm(egm_count,:)=value.map_signal;
                userdata.electric.elctrodeNames_uni(egm_count,:)=[value.electrode_number1, value.electrode_number2];
                userdata.electric.egmRef(egm_count,:)=value.ref_signal;
                userdata.electric.annotations.woi(egm_count,:)=[value.woi_start_index, value.woi_end_index]-value.ref_annotation_index;
                userdata.electric.annotations.referenceAnnot(egm_count,1)=value.ref_annotation_index;
                userdata.electric.annotations.mapAnnot(egm_count,1)=value.map_annotation_index;
                userdata.electric.egmSurfX(egm_count,:)=value.projected; %%This is the one on the surface
                
                
                %added to get index infomation for translating egm to shell
                userdata.electric.include(egm_count,1)=value.valid_for_map; %not discarded and has ref/map annotations so not necessairly just the opposite of valid_for_map/inlcude
                %userdata.electric.clipped(egm_count,1)=value.clipped;
                userdata.electric.discarded(egm_count,1)=value.discarded; %discarded means manual deletion of the point
                
            end
end

%ablation data
for i=1:numel(tempfolder_abl)
    fname=tempfolder_abl(i).name;
    [a,b,c]=fileparts(fname);
    
    
    if strcmp(c,'.json') == 1 && contains(b,'ablation')
                abl_count=abl_count+1;
                x=kodex_egm_import([abl_folder,filesep(),fname]);
                value = jsondecode(x);
                
                %% POPULATE ABLATION DATA HERE (rfindex)
                userdata.rfindex.tag.X(abl_count,:)=value.position;
                userdata.rfindex.tag.avgForce(abl_count,1)=value.average_cf;
                userdata.rfindex.tag.maxTemp(abl_count,1)=max(value.rf_temp_vec);
                userdata.rfindex.tag.maxPower(abl_count,1)=max(value.rf_power_vec);
                userdata.rfindex.tag.impedance.baseImp(abl_count,1)=value.first_imp;
                userdata.rfindex.tag.impedance.impDrop(abl_count,1)=value.first_imp-value.min_imp;
                userdata.rfindex.tag.fti(abl_count,1)=value.cf_time_integral_vec(end); %%this is an assumption that each point is the total 'so far': need to check 
                userdata.rfindex.tag.index.name(abl_count,1)=value.transmural_prediction; %is this correct?
                userdata.rfindex.tag.index.value(abl_count,1)=value.transmural_prediction;
                
                %time stuff
                t1=split(value.start_timestamp_str,'_');
                t1=split(t1{2},'-');
                t1_day = str2num(t1{1}); t1_hr = str2num(t1{2}); t1_sec=str2num(t1{3});
                t2=split(value.timestamp_finished,'_');
                t2=split(t2{2},'-');
                t2_day = str2num(t1{1}); t2_hr = str2num(t2{2}); t2_sec=str2num(t2{3});
                
                hr=(t2_hr-t1_hr)*60; sec=t2_sec-t1_sec;
                time=hr+sec;
                
                userdata.rfindex.tag.time(abl_count,1)=time;
    end
    
      %% POPULATE ABLATION DATA HERE (rf all data concentrated) - Need to work out
end

%landmark data
for i=1:numel(tempfolder_lm)
    fname=tempfolder_lm(i).name;
    [a,b,c]=fileparts(fname);
    if strcmp(c,'.json') == 1 && contains(b,'LM')
                lm_count=lm_count+1;
                x=kodex_egm_import([lm_folder,filesep(),fname]);
                value = jsondecode(x);
                %% POPULATE LANDMARK DATA HERE (added to electrogram data)
                
                userdata.electric.tags{egm_count+lm_count,1}=value.landmark_type;
                userdata.electric.names{egm_count+lm_count,1}=value.landmark_sub_type;
                userdata.electric.egmX(egm_count+lm_count,:)=value.position; %%is raw position
                userdata.electric.indexmesh(egm_count,1)=NaN;
                userdata.electric.voltages.bipolar(egm_count+lm_count,1)=NaN;
                userdata.electric.LATs(egm_count+lm_count,1)=-8000;
                userdata.electric.egm(egm_count+lm_count,:)=NaN(1,size(userdata.electric.egm,2));
                userdata.electric.elctrodeNames_uni(egm_count+lm_count,:)=NaN(1,size(userdata.electric.electrodeNames_uni,2));
                userdata.electric.egmRef(egm_count+lm_count,:)=NaN(1,size(userdata.electric.egmRef,2));
                userdata.electric.annotations.woi(egm_count+lm_count,:)=NaN(1,2);
                userdata.electric.annotations.referenceAnnot(egm_count+lm_count,1)=0;
                userdata.electric.annotations.mapAnnot(egm_count+lm_count,1)=0;
                userdata.electric.egmSurfX(egm_count+lm_count,:)=NaN(1,3); %%This is the one on the surface
                
                
                %added to get index infomation for translating egm to shell
                userdata.electric.include(egm_count+lm_count,1)=0;
                %userdata.electric.clipped(egm_count+lm_count,1)=NaN;
                userdata.electric.discarded(egm_count+lm_count,1)=1; %discarded will be opposite of include, can proabably remove
    end
end

                
end


if exist('study_dir_new') == 1
for i=1:numel(tempfolder);
    fname=tempfolder(i).name;
    [a,b,c]=fileparts(fname);
    if strcmp(c,'.json') == 1 && contains(b,chamber)
                egm_count=egm_count+1;
                x=kodex_egm_import([egm_folder,filesep(),fname]);
                value = jsondecode(x);
                
                index(egm_count)=value.index_on_mesh;
                
                userdata.electric.indexmesh(egm_count,1)=value.index_on_mesh;
                
                %electric data
                userdata.electric.tags{egm_count,1}=value.tag;
                userdata.electric.names{egm_count,1}=value.current_point;
                userdata.electric.egmX(egm_count,:)=value.position; %%is raw position
                userdata.electric.voltages.bipolar(egm_count,1)=Volts(value.current_point);
                userdata.electric.LATs(egm_count,1)=LATs(value.current_point);
                userdata.electric.egm(egm_count,:)=value.map_signal;
                userdata.electric.elctrodeNames_uni(egm_count,:)=[value.electrode_number1, value.electrode_number2];
                userdata.electric.egmRef(egm_count,:)=value.ref_signal;
                userdata.electric.annotations.woi(egm_count,:)=[value.woi_start_index, value.woi_end_index]-value.ref_annotation_index;
                userdata.electric.annotations.referenceAnnot(egm_count,1)=value.ref_annotation_index;
                userdata.electric.annotations.mapAnnot(egm_count,1)=value.map_annotation_index;
                userdata.electric.egmSurfX(egm_count,:)=value.projected; %%This is the one on the surface
                
                
                %added to get index infomation for translating egm to shell 
                userdata.electric.include(egm_count,1)=value.valid_for_map;
                userdata.electric.clipped(egm_count,1)=value.clipped;
                userdata.electric.discarded(egm_count,1)=value.discarded;
                
            end
end
        
%ablation data
for i=1:numel(tempfolder_abl)
    fname=tempfolder_abl(i).name;
    [a,b,c]=fileparts(fname);
    
    
    if strcmp(c,'.json') == 1 && contains(b,'ablation')
                abl_count=abl_count+1;
                x=kodex_egm_import([abl_folder,filesep(),fname]);
                value = jsondecode(x);
                
                %% POPULATE ABLATION DATA HERE (rfindex)
                userdata.rfindex.tag.X(abl_count,:)=value.position;
                userdata.rfindex.tag.avgForce(abl_count,1)=value.average_cf;
                userdata.rfindex.tag.maxTemp(abl_count,1)=max(value.rf_temp_vec);
                userdata.rfindex.tag.maxPower(abl_count,1)=max(value.rf_power_vec);
                userdata.rfindex.tag.impedance.baseImp(abl_count,1)=value.first_imp;
                userdata.rfindex.tag.impedance.impDrop(abl_count,1)=value.first_imp-value.min_imp;
                userdata.rfindex.tag.fti(abl_count,1)=value.cf_time_integral_vec(end); %%this is an assumption that each point is the total 'so far': need to check 
                userdata.rfindex.tag.index.name(abl_count,1)=value.transmural_prediction; %is this correct?
                userdata.rfindex.tag.index.value(abl_count,1)=value.transmural_prediction;
                
                %time stuff
                t1=split(value.start_timestamp_str,'_');
                t1=split(t1{2},'-');
                t1_day = str2num(t1{1}); t1_hr = str2num(t1{2}); t1_sec=str2num(t1{3});
                t2=split(value.timestamp_finished,'_');
                t2=split(t2{2},'-');
                t2_day = str2num(t1{1}); t2_hr = str2num(t2{2}); t2_sec=str2num(t2{3});
                
                hr=(t2_hr-t1_hr)*60; sec=t2_sec-t1_sec;
                time=hr+sec;
                
                userdata.rfindex.tag.time(abl_count,1)=time;
    end
    
      %% POPULATE ABLATION DATA HERE (rf all data concentrated) - Need to work out
end

%landmark data
for i=1:numel(tempfolder_lm)
    fname=tempfolder_lm(i).name;
    [a,b,c]=fileparts(fname);
    if strcmp(c,'.json') == 1 && contains(b,'LM')
                lm_count=lm_count+1;
                x=kodex_egm_import([lm_folder,filesep(),fname]);
                value = jsondecode(x);
                %assignin('base','LMvalue',value)
                %% POPULATE LANDMARK DATA HERE (added to electrogram data)
                
                userdata.electric.tags{egm_count+lm_count,1}=value.landmark_type;
                userdata.electric.names{egm_count+lm_count,1}=value.landmark_sub_type;
                userdata.electric.egmX(egm_count+lm_count,:)=value.position; %%is raw position
                userdata.electric.indexmesh(egm_count,1)=NaN;
                userdata.electric.voltages.bipolar(egm_count+lm_count,1)=NaN;
                userdata.electric.LATs(egm_count+lm_count,1)=-8000;
                userdata.electric.egm(egm_count+lm_count,:)=NaN(1,size(userdata.electric.egm,2));
                userdata.electric.elctrodeNames_uni(egm_count+lm_count,:)=NaN(1,size(userdata.electric.elctrodeNames_uni(1,:),2));
                userdata.electric.egmRef(egm_count+lm_count,:)=NaN(1,size(userdata.electric.egmRef,2));
                userdata.electric.annotations.woi(egm_count+lm_count,:)=NaN(1,2);
                userdata.electric.annotations.referenceAnnot(egm_count+lm_count,1)=0;
                userdata.electric.annotations.mapAnnot(egm_count+lm_count,1)=0;
                userdata.electric.egmSurfX(egm_count+lm_count,:)=NaN(1,3); %%This is the one on the surface
                
                
                %added to get index infomation for translating egm to shell
                userdata.electric.include(egm_count+lm_count,1)=0;
                %userdata.electric.clipped(egm_count+lm_count,1)=0;
                userdata.electric.discarded(egm_count+lm_count,1)=1; %discarded will be opposite of include, can proabably remove
    end
end


end


userdata.electric.numberofegms=egm_count;

if egm_count == 0 %This means no json file, so got to work just off csv file
    for i=1:size(csv_pos,1) %go through all positions
        egm_point=csv_pos(i,:);
        D=pdist2(egm_point,obj.v);
        [~,vertex]=min(D);
        index(i)=vertex;
        egm_vertex(i,:)=obj.v(vertex,:);
        userdata.electric.egmX(i,:)=egm_point;
        userdata.electric.egmSurfX(i,:)=egm_vertex(i,:);
        userdata.electric.voltages.bipolar(i,1)=Volts(i);
        userdata.electric.annotations.mapAnnot(i,1)=LATs(i);
        
    end
    userdata.electric.annotations.referenceAnnot(:,1)=zeros(size(userdata.electric.annotations.mapAnnot));
    
end

%% Get surface data
TR = TriRep(obj.f.vt, obj.v(:,1), obj.v(:,2), obj.v(:,3));
userdata.surface.triRep=TR;
%is vertex at rim (all 0 for these examples, but need to see if this will always be the case)
for i=1:size(obj.v,1)
    userdata.surface.isVertexAtRim(i)=false;
end

number_of_vertices = size(userdata.surface.triRep.X,1);

%% calcualte LATs and BiPs (mean of all values with same index on mesh)
% This is the method used by Kodex, possibly could be changed/updated?
userdata.surface.act_bip=NaN(number_of_vertices,2);
%first find 'extreme' indecies for exclusion from 'mean' measurements
[max_voltage,ind_temp]=max(Volts); ind_max_voltage=index(ind_temp);
%userdata.electric.voltages.bipolar(ind_max_voltage,1)=max_voltage;
userdata.surface.act_bip(ind_max_voltage,2)=max_voltage;

[min_voltage,ind_temp]=min(Volts); ind_min_voltage=index(ind_temp);
%userdata.electric.voltages.bipolar(ind_min_voltage,1)=min_voltage;
userdata.surface.act_bip(ind_min_voltage,2)=min_voltage;

[max_LAT,ind_temp]=max(LATs); ind_max_LAT=index(ind_temp);
userdata.surface.act_bip(ind_max_LAT,1)=max_LAT;

[min_LAT ind_temp]=min(LATs); ind_min_LAT=index(ind_temp);
userdata.surface.act_bip(ind_min_LAT,1)=min_LAT;

for i=1:numel(index);
    if i ~=  ind_max_voltage && i ~=  ind_min_voltage
        volts_index=Volts(index==i);
        if isempty(volts_index) == 0
            userdata.surface.act_bip(i,2)=mean(volts_index);
        end
    end
    
    if i ~=  ind_max_LAT && i ~=  ind_min_LAT
        lats_index=LATs(index==i);
        if isempty(lats_index) == 0
            userdata.surface.act_bip(i,1)=mean(lats_index);
        end
    end
end





%% Save data
if save_option == true
uisave('userdata')
end

savefilename
if isempty('savefilename') == 1
save(savefilename,'userdata')
end



