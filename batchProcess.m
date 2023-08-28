% ----------------------------------------------------------------------- %
% OPENEP/batchProcess is a template script which can be used to perform
% multiple functions on OpenEP datasets
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%                             Configuration
% Set up directories and load filenames
working_dir = '/Users/Steven/Desktop/openep_working_dir';
% ----------------------------------------------------------------------- %


% Get a list of files
disp('Getting list of filenames')
allFiles = nameFiles(working_dir, 'showhiddenfiles', false, 'extension', 'mat');

lat.carto = [];
lat.openep = [];

% Iterate through each file and perform some OpenEP functions
for i = 1:numel(allFiles)

    
   disp(['working on file ... ' allFiles{i}])
   load([working_dir filesep() allFiles{i}]);
   
   % Geometry
   openep_num_pts(i,1) = getNumPts(userdata);
   openep_chamber_area(i,1) = getArea(userdata);
   openep_chamber_area_closed(i,1) = getArea(userdata, 'method', 'fill');
   openep_chamber_volume(i,1) = getVolume(userdata);
   
   % Local activation time
   [~, earliestDirc(i,1:3)] = getEarliestActivationSite(userdata, 'method', 'ptbased');
   [~, earliestPrct(i,1:3)] = getEarliestActivationSite(userdata, 'method', 'ptbasedprct');
   tat.ptbased(i,1) = getTotalActivationTime(userdata, 'method', 'ptbased');
   tat.ptbasedprct(i,1) = getTotalActivationTime(userdata, 'method', 'ptbasedprct');
   tat.clinmap(i,1) = getTotalActivationTime(userdata, 'method', 'clinmap');
   tat.clinmapprct(i,1) = getTotalActivationTime(userdata, 'method', 'clinmapprct');
   tat.openepmap(i,1) = getTotalActivationTime(userdata, 'method', 'openepmap');
   tat.openepmapprct(i,1) = getTotalActivationTime(userdata, 'method', 'openepmapprct');
   lat_carto_temp = userdata.surface.act_bip(:,1);
   lat_openep_tmp = generateInterpData(userdata, 'lat-map');
   lat.carto = [lat.carto; lat_carto_temp];
   lat.openep = [lat.openep; lat_openep_tmp];
   
   % Voltage
   meanVoltage.Carto(i,1) = getMeanVoltage(userdata, 'method', 'map');
   meanVoltage.OpenEP(i,1) = getMeanVoltage(userdata, 'method', 'egm');
   lva.Carto(i,1) = getLowVoltageArea(userdata, 'method', 'map');
   lva.OpenEP(i,1) = getLowVoltageArea(userdata, 'method', 'egm');
   
   % Conduction Velocity
end