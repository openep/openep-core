% kodex_to_openEP
% converter for kodex data to poptulate a data strucutre that can be read by openEP
%COS 22/12 - Currently works for kodex data of strcutre Folder -> studies -> 'LAT' 
%'LAT' folder then contains csv points file (lat_points), mesh file (.obj) and (optional) .json files with raw electrogram data 

clear all

%% Upload 
study_name1=uigetdir;
study_name=[study_name1];
xx=dir(study_name);
map_count=0;
for j=3:numel(xx)
if xx(j).isdir == 1
   map_count=map_count+1;
   studies{map_count}=xx(j).name;
end
end

[indx,tf] = listdlg('Liststring',studies);
% indx=1;
chamber=char(studies(indx));
egm_folder=[study_name,'\',chamber,'\LAT'];
xx=dir(egm_folder);




mesh_folder=[study_name1,'\',chamber,'\','LAT'];
files=dir(mesh_folder);

Fs=500; %work out way to read frequency
sampling_rate=1000/Fs; %in ms

%% find object and report data files (as measured by kodex)
for n=3:numel(files)
    fname=files(n).name;
    [a,b,c]=fileparts(fname);
    if strcmp('.obj',c) == 1
    object_file=[mesh_folder,'\',b,c];
    end
    if strcmp('.csv',c) == 1 && strcmp('lat_points',b) == 1
    csv_file=[mesh_folder,'\',b,c];
    end
    if strcmp('.csv',c) == 1 && strcmp('landmark_points',b) == 1
    csv_landmark_file=[mesh_folder,'\',b,c];
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
%% this is to fill in from 'raw' electrogram data if avalable
for i=1:numel(xx);
   fname=xx(i).name;
   [a,b,c]=fileparts(fname);
   if strcmp(c,'.json') == 1
      
       
      egm_count=egm_count+1;
      x=kodex_egm_import([egm_folder,'\',fname]);
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
    userdata.surface.act_bip(i,2)=mean(volts_index)
    end
    end
    
    if i ~=  ind_max_LAT && i ~=  ind_min_LAT
    lats_index=LATs(index==i);
    if isempty(lats_index) == 0
    userdata.surface.act_bip(i,1)=mean(lats_index)
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

uisave('userdata')

clearvars -except userdata 


