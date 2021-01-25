function cemrg2carto(infile, outfile, vtktitle)
% CEMRG2CARTO Converts a VTK for loading into Carto
% Usage:
%   op = cemrg2carto(infile, outfile)
%
% Where:
%   infile - a VTK file created by CEMRG
%   outfile - a VTK file formatted for Carto
%
% CEMRG2CARTO Converts a VTK file created by CEMRG into a VTK file that can
% be loaded into Carto. The two main functions are to convert cell data to
% point data and to add the correct header based on patient name and ID.
%
% Author: Steven Williams (2016)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% 
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

if strcmpi(infile, 'openfile')
    [filename, pathname] = uigetfile('*.vtk', 'Select VTK file to open ...');
    infile = [pathname filename];
end
if infile==0
    disp('Operation cancelled');
    return
end

h = waitbar(0, 'Load VTK file');
hVtk = VTKReader(infile);

waitbar(.2, h, 'Read VTK data');
hVtk.readAllData();

waitbar(.5, h, 'Get anatomy from VTK file');
tr = hVtk.getTriRep();

waitbar(.8, h, 'Convert cell data to point data');
pointdata = trFaceToVertData(tr, hVtk.CellData{1});
close(h);

if exist('vtktitle')==0
    answer = inputdlg({'firstname', 'surname', 'ID'});
    vtktitle = ['PatientData ' answer{1} ' ' answer{2} ' ' answer{3}];
end

if strcmpi(outfile, 'openfile')
    [filename, pathname] = uiputfile('*.vtk', 'Select VTK file to save ...');
    outfile = [pathname filename];
end

writeTriRep2VTK(tr, pointdata, 'title', vtktitle, 'outputfile', outfile);