% ----------------------------------------------------------------------- %
% OPENEP/openep-demo is a demo script to test the installation
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%                             Configuration
% Configuration info goes here
% ----------------------------------------------------------------------- %

disp('Loading OpenEP Dataset 1 ...');
load openep_dataset_1

disp('Create a figure ...');
figure
for i = 1:4
    hAx(i) = subplot(2,2,i);
end

disp('Conduction velocity mapping ...')
cvdata = getConductionVelocity(userdata);
cvdata(cvdata>2) = NaN;
axes(hAx(1));
drawMap(userdata, 'type', 'cv', 'orientation', 'ap', 'data', cvdata);

disp('Finding anatomical structures ...')
axes(hAx(3));
FF = getAnatomicalStructures(userdata, 'plot', true);

disp('Show a voltage map ...')
axes(hAx(2));
drawMap(userdata, 'type', 'bip', 'coloraxis', [0 2]);

disp('Voltage histogram analsysis')
axes(hAx(4));
thresh = [ 0.01 0.5; 0.51 1; 1.01 1.50; 1.51 2; 2.01 2.50 ];
area = voltageHistogramAnalysis(userdata, 'plot', 'true', 'threshold', thresh);

close(gcf);