function op = writeTriRep2VTK(tr, pointdata, varargin)
% WRITETRIREP2VTK Converts a triRep to a VTK.
% Usage:
%   op = writeTriRep2VTK(tr, pointdata)
%   op = writeTriRep2VTK(tr, pointdata, 'outputfile', '<pathttovtk>')
% Where:
%   tr - is the triRep
%   pointdata - is the point data to be written to the VTK, or []
%
% WRITETRIREP2VTK Converts a triRep to a VTK file. Scalar point data is
% specified. Param-Value pairs:
%   'outputfile' - 'openfile' or absolute path
%   'scalarnames' - cell array of names. Number of elements given by
%           size(pointdata,2);
%   'type' - 'ascii' or 'binary' 
%   'celldata' - optional cell data to write to the VTK file
%   'title' - optional title to write to the VTK file
%
% Author: Steven Williams (2013)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% test code
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

outputfile = 'openfile';
scalarnames = [];
type = 'ascii';
celldata = [];
title = '';

if nargin==1
    pointdata = [];
end

if nargin>2
    for i = 1:2:nargin-2 %2 static arguments
        switch(varargin{i})
            case 'outputfile'
                outputfile = varargin{i+1};
            case 'scalarnames'
                scalarnames = varargin{i+1};
            case 'type'
                type = varargin{i+1};
            case 'celldata'
                celldata = varargin{i+1};
            case 'title'
                title = varargin{i+1};
        end
    end
end

if strcmpi(outputfile, 'openfile')
    [filename, pathname] = uiputfile('*.vtk', 'Pick a VTK file');
    outputfile = [pathname filename];
end
if isempty(scalarnames) && ~isempty(pointdata)
   scalarnames = cell(size(pointdata,2),1);
   for i = 1:numel(scalarnames)
       scalarnames{i} = ['scalars' num2str(i)];
   end
end

% open the file for writing
fid = fopen(outputfile, 'w');

% write the header
fprintf(fid, '%s\n', '# vtk DataFile Version 3.0');
if isempty(title)
    fprintf(fid, '%s\n', 'vtk output');
else
    fprintf(fid, '%\n', title);
end
switch type
    case 'ascii'
        fprintf(fid, '%s\n', 'ASCII');
    case 'binary'
        fprintf(fid, '%s\n', 'BINARY');
end

% write the polygon data
fprintf(fid, '%s\n', 'DATASET POLYDATA');
numpts = size(tr.X, 1);
fprintf(fid, '%s\t', 'POINTS');
fprintf(fid, '%d\t', numpts);
fprintf(fid, '%s\n', 'float');
switch type
    case 'ascii'
        fprintf(fid, '%f %f %f\n', tr.X');
    case 'binary'
        fwrite(fid, [reshape(tr.X(:,1),1,numpts);  reshape(tr.X(:,2),1,numpts); reshape(tr.X(:,3),1,numpts)], 'float', 'b'); %BINARY
end
fprintf(fid, '\n');

polygons = zeros(size(tr.Triangulation,1),4);
polygons(:,1) = 3;
polygons(:,2:4) = tr.Triangulation - 1;

numtri = size(tr.Triangulation, 1);
fprintf(fid, '%s\t', 'POLYGONS');
fprintf(fid, '%d\t', numtri);
fprintf(fid, '%d\n', numtri*4);
switch type
    case 'ascii'
        fprintf(fid, '%d %d %d %d\n', polygons');
    case 'binary'
        fwrite(fid, [reshape(polygons(:,1),1,numtri); reshape(polygons(:,2),1,numtri); reshape(polygons(:,3),1,numtri); reshape(polygons(:,4),1,numtri)], 'int', 'b'); %BINARY
end

if ~isempty(pointdata)
    fprintf(fid, '\n%s\t', 'POINT_DATA');
    numdatapts = size(pointdata,1);
    fprintf(fid, '%d\n', size(pointdata,1));
    if size(pointdata,2)==1
        i=1;
        fprintf(fid, '%s\n', ['SCALARS scalars float']);
        fprintf(fid, '%s\n', 'LOOKUP_TABLE default');
        switch type
            case 'ascii'
                pointdata(isnan(pointdata)) = -999;
                fprintf(fid, '%f ', pointdata(:,i));
            case 'binary'
                fwrite(fid, reshape(pointdata(:,i),1,numdatapts), 'float', 'b'); %BINARY was 'float'
        end
        fprintf(fid, '\n');
    else
        
        for i = 1:size(pointdata,2)
            fprintf(fid, '%s\n', ['SCALARS ' scalarnames{i} ' float']);
            fprintf(fid, '%s\n', 'LOOKUP_TABLE default');
            switch type
                case 'ascii'
                    pointdata(isnan(pointdata)) = -999;
                    fprintf(fid, '%f ', pointdata(:,i));
                case 'binary'
                    fwrite(fid, reshape(pointdata(:,i),1,numdatapts), 'float', 'b'); %BINARY was 'float'
            end
            fprintf(fid, '\n');
        end
    end
end

if ~isempty(celldata)
    disp('cell data written')
    fprintf(fid, '\n%s\t', 'CELL_DATA');
    numdatapts = size(celldata,1);
    fprintf(fid, '%d\n', size(celldata,1));
    for i = 1:size(celldata,2)
        %fprintf(fid, '%s\n', ['SCALARS ' 'CellData' ' float']);  
        fprintf(fid, '%s\n', ['SCALARS ' scalarnames{i} ' float']);
        fprintf(fid, '%s\n', 'LOOKUP_TABLE default');
        switch type
            case 'ascii'
                celldata(isnan(celldata)) = -999;
                fprintf(fid, '%f ', celldata(:,i));
            case 'binary'
                fwrite(fid, reshape(celldata(:,i),1,numdatapts), 'float', 'b'); %BINARY was 'float'
        end
        fprintf(fid, '\n');
    end
end

op = true;
end