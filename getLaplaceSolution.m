function selpath = getLaplaceSolution (tvIndices, ivcIndices, svcIndices, path,userdata,rootMeshFile)
% GETULAPLACESOLUTION Returns the solved laplace solution for UAC
% calculation. It calls openCarp to solve the Laplace for the given
% geometries and Paths provided.
%
% Usage:
%   getLaplaceSolution (tvIndices, ivcIndices, svcIndices, path)
% Where:
%   userdata   - an OpenEP data structure
%
% GETUACBOUNDARYCONDITIONS accepts the following parameter-value pairs
%   'plot'     {false}|true
%
% GETUACBOUNDARYCONDITIONS identifies the boundary conditions for UAC calculation.
%
% Author: Ali Gharaviri (2022) (Copyright)
% 
%
% Modifications -
%
% See also FREEBOUNDARYPOINTS
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------


%%%%%%%%% writing files required to run laplace operation on openCarp

selpath = uigetdir(pwd,'select a directory to write Carp files');
cd (selpath)
%%% 1- Write elem and pts files
writeOpenEP2CarpMesh(userdata, rootMeshFile)

%%%% 2- write Lon file
fidend = fopen([rootMeshFile,'.lon'],'w');
lonData = [ones(size(userdata.surface.triRep.Triangulation,1),1) zeros(size(userdata.surface.triRep.Triangulation,1),2)];
fprintf(fidend,'%d %d %d',lonData(1,1:3));
for numLonElem = 2 : size(lonData,1)
    fprintf(fidend,'\n%d %d %d',lonData(numLonElem,1:3));
end

%%%%% 3- Defining borders to calculate Laplace
% 3-1- Border 1: TCV valve
fidend = fopen(['TCV','.vtx'],'w');

fprintf(fidend,'%d',size(tvIndices,1));
fprintf(fidend,'\n%s','intra');
for numIndTV=1:size(tvIndices,1)
    fprintf(fidend,'\n%d',tvIndices(numIndTV,1)-1);
end
fclose(fidend);


%%%%% 3-2 Boarder 2: ICV-SCV
fidend = fopen(['ICV_SCV','.vtx'],'w');
fprintf(fidend,'%d',size(path{1,1},1));
fprintf(fidend,'\n%s','intra');
for numIndPath=1:size(path{1,1},1)
    fprintf(fidend,'\n%d',path{1,1}{numIndPath,1}.id-1); 
end
fclose(fidend);
%%%%% 3-3 Boarder 3: ICV
fidend = fopen(['ICV','.vtx'],'w');
fprintf(fidend,'%d',size(ivcIndices,1));
fprintf(fidend,'\n%s','intra');
for numIndIVC=1:size(ivcIndices,1)
    fprintf(fidend,'\n%d',ivcIndices(numIndIVC,1)-1);
end
fclose(fidend);

%%%%% 3-4 Border 4: SCV
fidend = fopen(['SCV','.vtx'],'w');
fprintf(fidend,'%d',size(svcIndices,1));
fprintf(fidend,'\n%s','intra');
for numIndSVC=1:size(svcIndices,1)
    fprintf(fidend,'\n%d',svcIndices(numIndSVC,1)-1);%%%%% Be careful just tested
end

fclose(fidend);

%%%%% 3 -5 Change .par files
selpath = uigetdir(pwd,'select a directory with sample .par files');
cd (selpath)

% 3-5-1: write par file for get laplace sollution for Alpha Coordinate
parFile = readlines('sample_1.par');
parFile (1) = sprintf('meshname = %s', rootMeshFile);
fidend = fopen(['RA_3','.par'],'w');
fprintf(fidend,'%s\n',parFile);
fclose(fidend)
% 3-5-1: write par file for get laplace sollution for Beta Coordinate
parFile = readlines('sample_2.par');
parFile (1) = sprintf('meshname = %s', rootMeshFile);
fidend = fopen(['RA_4','.par'],'w');
fprintf(fidend,'%s\n',parFile);
fclose(fidend)

%%%% Change directory where you want to load input files and save openCarp outcomes
selpath = uigetdir(pwd,'select a directory to write .par files');

movefile ('RA_3.par', selpath)
movefile ('RA_4.par', selpath)
%movefile ('*.vtx', selpath)
%movefile ('*.elem', selpath)
%movefile ('*.lon', selpath)
%movefile ('*.pts', selpath)

cd (selpath)



%%%%% Calling openCarp from Matlab and run the Laplace operaion

!openCARP +F ./RA_3.par 
cd ./OUTPUT_DIR/
phieFileName = sprintf('%s%s',rootMeshFile,'_phie_alpha.igb');
phie_i_FileName = sprintf('%s%s',rootMeshFile,'_phie_i_alpha.igb');
 movefile ('phie.igb', phieFileName);
 movefile ('phie_i.igb', phie_i_FileName);

cd ..

!openCARP +F ./RA_4.par 
cd ./OUTPUT_DIR/
betaFileName = sprintf('%s%s',rootMeshFile,'_phie_beta.igb');
beta_i_FileName = sprintf('%s%s',rootMeshFile,'_phie_i_beta.igb');
movefile ('phie.igb', betaFileName)
movefile ('phie_i.igb', betaFileName)


%end