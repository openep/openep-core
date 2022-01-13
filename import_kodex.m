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

%% Upload

if nargin == 0
    study_dir_overall=uigetdir;
else
    study_dir_overall=varargin{1}
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

map_count=0;
for j=3:numel(tempfolder)
    if tempfolder(j).isdir == 1
        map_count=map_count+1;
        studies{map_count}=tempfolder(j).name;
    end
end

%check format and if study directory needs to be changed to 'full' format
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
    tempfolder=dir(egm_folder);
    mesh_folder=[study_dir_overall,filesep(),chamber,filesep(),'LAT'];
    files=dir(mesh_folder);
elseif exist('study_dir_new') == 1
    if isempty(chamber) == 1
        [indx,tf] = listdlg('Liststring',studies);
        chamber=char(studies(indx))
    end
    egm_folder=[study_dir_new,filesep(),chamber,filesep(),'LAT'];
    tempfolder=dir(egm_folder);
    mesh_folder=[study_dir_overall,filesep(),'Report',filesep(),chamber];
    files=dir(mesh_folder);
end


Fs=500; %work out way to read frequency
sampling_rate=1000/Fs; %in ms

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
datas=split(table2array(Values),',');
LATs=str2double(datas(:,5));
LATs=LATs-min(LATs);
Volts=str2double(datas(:,6));
csv_pos=str2double(datas(:,2:4));
obj = readObj(object_file);


T1=readtable(csv_landmark_file);
W1 = width(T1);
Landmarks=T1(:,W1);
datas2=split(table2array(Landmarks),',');
Lpoints=str2double(datas2(:,2:4));

%% Get 'raw' electirc data
egm_count=0;
userdata = openep_createuserdata();
%% this is to fill in from 'raw' electrogram data if avalable (think of better way)

if exist('study_dir_new') == 0
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
                userdata.electric.names(egm_count,1)=value.current_point;
                userdata.electric.egmX(egm_count,:)=value.position; %%is raw position
                userdata.electric.voltages.bipolar(egm_count,1)=Volts(value.current_point);
                userdata.electric.LATs(egm_count,1)=LATs(value.current_point);
                userdata.electric.egm(egm_count,:)=value.map_signal;
                userdata.electric.elctrodeNames_uni(egm_count,:)=[value.electrode_number1, value.electrode_number2];
                userdata.electric.egmRef(egm_count,:)=value.ref_signal;
                userdata.electric.annotations.woi(egm_count,:)=[value.woi_start_index, value.woi_end_index];
                userdata.electric.annotations.referenceAnnot(egm_count,1)=value.ref_annotation_index;
                userdata.electric.annotations.mapAnnot(egm_count,1)=value.map_annotation_index;
                userdata.electric.egmSurfX(egm_count,:)=value.projected; %%This is the one on the surface
                
                
                %added to get index infomation for translating egm to shell (PROBABLY WRONG)
                userdata.electric.valid_for_map(egm_count,1)=value.valid_for_map;
                userdata.electric.clipped(egm_count,1)=value.clipped;
                userdata.electric.discarded(egm_count,1)=value.discarded;
                
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
                userdata.electric.names(egm_count,1)=value.current_point;
                userdata.electric.egmX(egm_count,:)=value.position; %%is raw position
                userdata.electric.voltages.bipolar(egm_count,1)=Volts(value.current_point);
                userdata.electric.LATs(egm_count,1)=LATs(value.current_point);
                userdata.electric.egm(egm_count,:)=value.map_signal;
                userdata.electric.elctrodeNames_uni(egm_count,:)=[value.electrode_number1, value.electrode_number2];
                userdata.electric.egmRef(egm_count,:)=value.ref_signal;
                userdata.electric.annotations.woi(egm_count,:)=[value.woi_start_index, value.woi_end_index];
                userdata.electric.annotations.referenceAnnot(egm_count,1)=value.ref_annotation_index;
                userdata.electric.annotations.mapAnnot(egm_count,1)=value.map_annotation_index;
                userdata.electric.egmSurfX(egm_count,:)=value.projected; %%This is the one on the surface
                
                
                %added to get index infomation for translating egm to shell (PROBABLY WRONG)
                userdata.electric.valid_for_map(egm_count,1)=value.valid_for_map;
                userdata.electric.clipped(egm_count,1)=value.clipped;
                userdata.electric.discarded(egm_count,1)=value.discarded;
                
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

%% calcualte LATs and BiPs (mean of all values with same index on mesh)

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

userdata.surface.act_bip(userdata.surface.act_bip==0)=NaN;


%% Get surface data
TR = TriRep(obj.f.vt, obj.v(:,1), obj.v(:,2), obj.v(:,3));
userdata.surface.triRep=TR;
%is vertex at rim (all 0 for these examples, but need to see if this will always be the case)
for i=1:size(obj.v,1)
    userdata.surface.isVertexAtRim(i)=false;
end


%% Save data
if save_option == true
uisave('userdata')
end

savefilename
if isempty('savefilename') == 1
save(savefilename,'userdata')
end



