% kodex_to_openEP
% converter for kodex data to poptulate a data strucutre that can be read by openEP
clear all
load('cartostruct.mat')
%% Upload 
study_name1=uigetdir;
%study_name1='2021-04-15_09-42-57_HR__2021-04-15_10-18-23'
study_name=[study_name1,'\dump_state\MainEntry\'];
xx=dir(study_name);
map_count=0;
for j=3:numel(xx)
if xx(j).isdir == 1
   map_count=map_count+1;
   studies{map_count}=xx(j).name;
end
end

% [indx,tf] = listdlg('Liststring',studies);
indx=1;
chamber=char(studies(indx));
egm_folder=[study_name,chamber,'\LAT'];
xx=dir(egm_folder);

mesh_folder=[study_name1,'\Report\',chamber];
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
datas=split(T.Var6,',');
LATs=str2double(datas(:,5));
%LATs=LATs-min(LATs);
Volts=str2double(datas(:,6));
csv_pos=str2double(datas(:,2:4));
obj = readObj(object_file);
userdata.surface.triRep=TriRep(obj.f.v,obj.v);

T1=readtable(csv_landmark_file);
datas2=split(T1.Var5,',');
Lpoints=str2double(datas2(:,2:4));

%% Get 'raw' electirc data
egm_count=0;
for i=1:numel(xx);
   fname=xx(i).name;
   [a,b,c]=fileparts(fname);
   if strcmp(c,'.json') == 1 && ~contains(b,chamber)==0
      
       
      egm_count=egm_count+1;
      x=kodex_egm_import([egm_folder,'\',fname]);
      value = jsondecode(x);
      
      
      %electric data
      userdata.electric.tags{egm_count,1}=value.tag;
      userdata.electric.names(egm_count,1)=value.current_point;
      userdata.electric.electrodeNames_bip{egm_count}='NaN'; %Bipolar electrode names?
      userdata.electric.egmX(egm_count,:)=value.position; %%assume this is raw position
      userdata.electric.egm(egm_count,:)=value.map_signal;
      userdata.electric.elctrodeNames_uni(egm_count,:)=[value.electrode_number1, value.electrode_number2];
      userdata.electric.egmUniX(egm_count,1)=NaN;
      userdata.electric.egmUni(egm_count,1)=NaN;
      userdata.electric.egmRef(egm_count,:)=value.ref_signal;
      userdata.electric.ecg(egm_count)=NaN;
      userdata.electric.annotations.woi(egm_count,:)=[value.woi_start_index, value.woi_end_index];
      userdata.electric.annotations.referenceAnnot(egm_count,1)=value.ref_annotation_index;
      userdata.electric.annotations.mapAnnot(egm_count,1)=value.map_annotation_index;
      userdata.electric.voltages.biploar(egm_count,1)=Volts(egm_count); 
      userdata.electric.voltages.uniploar(egm_count,1)=NaN; 
      userdata.electric.impedances.time(egm_count,1)=NaN;
      userdata.electric.impedances.value(egm_count,1)=NaN;
      %userdata.electric.egmSurfX(egm_count,:)=value.projected; %%This is the one on the surface
      userdata.electric.barDirection(egm_count,:)=[NaN,NaN,NaN]; 
      
      %added to get index infomation for translating egm to shell (PROBABLY WRONG)
      userdata.electric.indexmesh(egm_count,1)=value.index_on_mesh;
      userdata.electric.egmSurfX(egm_count,:)=obj.v(value.index_on_mesh,:);
         end
end
userdata.electric.numberofegms=egm_count;

%Get egm poistion data projected on to surface (find nearest point) PROBABLY WRONG
userdata.electric.egmX=csv_pos;
egm_vertex = zeros(size(csv_pos,1),3);

for i=1:size(csv_pos,1)
    egm_point=csv_pos(i,:);
    D=pdist2(egm_point,obj.v);
    [~,vertex]=min(D);
    egm_vertex(i,:)=obj.v(vertex,:);
end
    
userdata.electric.egmSurfX=egm_vertex;




%% Get surface data

userdata.surface.act_bip=NaN(size(obj.v,1),2);
userdata.surface.uni_imp_frc=NaN(size(obj.v,1),3);

%is vertex at rim (all 0 for these examples, but need to see if this will always be the case)
for i=1:size(obj.v,1)
    userdata.surface.isVertexAtRim(i)=false;
end

%match egm data to vertex
for i=1:egm_count;
    xx=userdata.electric.indexmesh(i);
    userdata.surface.act_bip(xx,1)=userdata.electric.annotations.referenceAnnot(i)-userdata.electric.annotations.mapAnnot(i); %this might need to be reversed as ref is after map sig
    userdata.surface.act_bip(xx,2)=Volts(i);
    userdata.surface.uni_imp_frc(xx,1)=userdata.electric.voltages.uniploar(i);
    userdata.surface.uni_imp_frc(xx,2)=userdata.electric.impedances.time(i);
    userdata.surface.uni_imp_frc(xx,3)=NaN;
end

%%
userdata.electric.annotations.referenceAnnot(egm_count,1)=value.ref_annotation_index;

% out_fac=1;
% %k=211;
% figure,
% hold on
% hSurf = trisurf(userdata.surface.triRep, 'edgecolor', 'none');
% axis equal vis3d
% set(hSurf, 'facecolor', [.5 .5 .5]);
% alpha(.3)
% set(gcf, 'color', 'white');
% axis off;
% hold on
% %plot3(Lpoints(:,1),Lpoints(:,2),Lpoints(:,3),'.r','MarkerSize',20)
% %  plot3(userdata.electric.egmSurfX(k,1),userdata.electric.egmSurfX(k,2),userdata.electric.egmSurfX(k,3),'.r','MarkerSize',20)
% % plot3(userdata.electric.egmX(k,1),userdata.electric.egmX(k,2),userdata.electric.egmX(k,3),'.b','MarkerSize',20)
% % plot3(csv_pos(k,1),csv_pos(k,2),csv_pos(l,3),'.g','MarkerSize',20)
% %plot3(userdata.electric.egmSurfX(:,1),userdata.electric.egmSurfX(:,2),userdata.electric.egmSurfX(:,3),'.r')
% %plot3(userdata.electric.egmX(:,1),userdata.electric.egmX(:,2),userdata.electric.egmX(:,3),'.b')
% plot3(csv_pos(:,1),csv_pos(:,2),csv_pos(:,3),'.g','MarkerSize',20)
% %plot3(egm_vertex(:,1),egm_vertex(:,2),egm_vertex(:,3),'.b')

%% Do some manual analyis of the electrograms
%kodex_analysis

