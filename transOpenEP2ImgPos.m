function [Y, iX, colourIndex, colourMap] = transOpenEP2ImgPos(userdata, openEp3DMesh, openEpUacMesh, img3DMesh, imgUacMesh, varargin)
% TRANSOPENEP2IMGPOS Translates electrode positions from OpenEP space to Imaging space
%
% Usage:
%   [Y, Ix] = transOpenEP2ImgPos(userdata, openEp3DMesh, openEpUacMesh, img3DMesh, imgUacMesh, varargin)
% Where:
%   userdata      - an OpenEP dataset
%   openEp3DMesh  - on OpenEP geometry, presented as a mesh structure,
%                   (see io_readCARPMesh.m)
%   openEpUacMesh - an OpenEP geometry converted into UACs (universal atrial
%                   co-ordinates), presented as a mesh structure
%   img3DMesh     - an imaging geometry, presented as a mesh structure
%   imgUacMesh    - an imaging geometry converted into UACs, presented as
%                   a mesh structure
%   Y             - 3D Cartesian co-ordinates of the electrode positions in imaging space
%   iX            - indexes into userdata.electric.X and identifies the
%                   electrodes which were translated:
%                       userdata.electric.X(iX,:) => Y
%   colourIndex   - an array of values indexing into colourMap, one value
%                   per co-ordinate triple inY
%   colourMap     - a 2D RGB colour map providing unique colours to
%                   electrodes based on their UACs. 
%
% TRANSOPENEP2IMGPOS accepts the following parameter-value pairs
%   'plot'         {false} | true
%           - specifies whether to draw a figure illustrating the process
%   'elecsamples'  {[]} | integer array
%           - integer array identifying the electrodes in
%             userdata.electric.X which will be drawn, if not empty
%   'eleccolors' {[]} | colorspec
%           - colors to draw the electrodes in the figures. Must satisfy
%           the condition size(eleccolors) == [numel(elecsamples),3]
%
% TRANSOPENEP2IMGPOS Uses the Universal Atrial Co-ordinates to convert
% electrode positions between OpenEP and Imaging co-ordinate systems. For
% example; this function can be used to create surface and 3D co-ordinates
% representing all of the recording electrode positions during a clinical
% case, relative to an imaging dataset such as MRI or CT.
%
% This translation is performed using the UACs to give the "surface-point-to
% -surface-point" translation between OpenEP and Imaging datasets. In order
% to determine the distance from the surface, i.e. how far internal to the
% surface the new points should be in the Imaging space, the following
% algorithm is used.
%   1. Lines are drawn between every vertex of the OpenEP mesh and the
%      barycentre of the OpenEP mesh.
%   2. For each electrode, the closest line is identified and the position
%      of the electrode is projected perpendicularly onto this line,
%      bisecting the line.
%   3. The ratio of the distance from the surface to the bisection (S1) to the
%      total length of the line (S2) is calculated.
%   4. Every line is uniquely identifiable by its surface vertex
%      attachment. The indices of these vertices is used to convert 3D
%      locations into 2D locations (openEp3DMesh -> openEpUacMesh)
%   5. Vertex indices in the MRI mesh are then identified by finding the
%      closest vertices in the imaging unfold mesh (openEpUacMesh -> imgUacMesh)
%   6. The barycentre of the imaging mesh is calculated and the lengths of
%      these lines are calculated (S3).
%   7. Finally, every electrode is translated the relevant distance along
%      the relevant surface-to-barycentre line, where this distance is
%      given by S3 * S1/S2
%
% Author: Steven Williams (2021)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
%  load('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Carto/Study_1_07_02_2017_19-42-07_1-Map.mat');
%  mOpenEp3D = io_readCarpMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Carto/Labelled');
%  mOpenEpUac = io_readCarpMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Carto/Labelled_Coords_2D_Rescaling_v3_C');
%  mImg3D = io_readCarpMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Model_16/Labelled');
%  mImgUac = io_readCarpMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Model_16/Labelled_Coords_2D_Rescaling_v3_C');
%  [Y, iX] = transOpenEP2ImgPos(userdata, mOpenEp3D, mOpenEpUac, mImg3D, mImgUac, 'plot', true, 'elecsamples', [101 529], 'eleccolors', [colorBrewer(1); colorBrewer(2)]);
%  [Y, iX, colourIndex, colourMap] = transOpenEP2ImgPos(userdata, mOpenEp3D, mOpenEpUac, mImg3D, mImgUac, 'plot', false, 'elecsamples', ':');
%  writeOpenEPElec2VTKImgPositions(Y, colourIndex, colourMap)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 5;
plot = false;
elecsamples = ':';
eleccolors = [];
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'plot'
                plot = varargin{i+1};
            case 'elecsamples'
                elecsamples = varargin{i+1};
            case 'eleccolors'
                eleccolors = varargin{i+1};
        end
    end
end

% Load the co-ordinates of the electrode positions and identify electrodes internal to the imaging mesh
X_orig = getElectrogramX(userdata);
tr = getTriangulationFromMeshStruct(openEp3DMesh, 'region', 11, 'scale', 'um', 'type', 'trirep');
iX = getMappingPointsWithinMesh(userdata, 'mask', tr);
X = X_orig(iX,:); % X is now all the points of interest but DOES NOT index into userdata.electric

% 1. Find the barycentre lines connecting vertices to C_ep
tempUserdata.surface.triRep = tr;
C_ep = getCentreOfMass(tempUserdata);
tr = getOpenEPSubsetFromImgMap(userdata, openEp3DMesh, 11);
vertices = tr.X;

% 2 & 3. For every electrode find the closest barycentre line (relevantVertex) and the  distance along this line from the OpenEP surface (distanceRatio)
for i = 1:length(X)
    distance = NaN(length(vertices),1);
    for j = 1:length(vertices)
        distance(j) = point_to_line_segment_distance(X(i,:), vertices(j,:), C_ep);
    end
    [~, iVertex4Electrode(i)] = nanmin(distance);
    lineSegmentLength(i) = lineLength([vertices(iVertex4Electrode(i),:); C_ep]);
    surf2ElectrodeDistance(i) = lineLength(    [vertices(iVertex4Electrode(i),:); X(i,:)]    );
    rV(i,1:3) = vertices(iVertex4Electrode(i),:); %rV for relevantVertex
    dR(i) = surf2ElectrodeDistance(i) / lineSegmentLength(i); %dR for distanceRatio
end

% 4. Convert 3D locations to 2D locations (openEp3DMesh -> openEpUacMesh)
iVertex_openEp3D = findclosestvertex(openEp3DMesh.Pts/1000, rV);
coordUac = openEpUacMesh.Pts(iVertex_openEp3D,:);

% 5. Find the closest vertex to this location in the MRI-UAC mesh; this defines the barycentre line of interest and sub
iVertex_mImgUac = findclosestvertex(imgUacMesh.Pts, coordUac);
nV = img3DMesh.Pts(iVertex_mImgUac,:); %nV for newvertex

% 6. Calculate the barycentre of the MRI-3D mesh, the MRI barycentre
% vectors and their lengths
tempUserdata.surface.triRep = getTriangulationFromMeshStruct(img3DMesh, 'type', 'trirep');
C_img = getCentreOfMass(tempUserdata);
baryVectors = repmat(C_img, size(nV,1), 1) - nV;
baryVectorLengths = vecnorm(baryVectors,2,2);

% 7. Translate UAC positions along the new vectors.
uV = baryVectors./baryVectorLengths; % uV for unitVector
scaledVectors = uV .* repmat(dR', 1, 3) .* baryVectorLengths;
Y = nV + scaledVectors;
X = X * 1000;

% work out the electrode color maps using a 2D color map
if strcmpi(elecsamples, ':')
    elecsamples = [1:1:size(X,1)];
end

R=[1 0; 
    1 0];
G=[1 1
    0 0];
B=[0 0
    0 1];
R = interp2(R,8);
G = interp2(G,8);
B = interp2(B,8);
I = 255*cat(3,R,G,B);
colourMap = repack2DColorMap(I);
for i = 1:numel(elecsamples)
    A = floor(1+coordUac(elecsamples(i),1)*255);
    B = floor(1+coordUac(elecsamples(i),2)*255);
    colors{i} = squeeze(I(A,B,:))'/256;
    colourIndex(i) = repack2DData(A,B,size(I,2));
end

if ~isempty(eleccolors)
   if size(eleccolors,1) ~= numel(elecsamples)
       error('OPENEP/TRANSOPENEP2IMGPOS: mismatch between the size of the values for eleccolors and elecsamples');
   else
       for i = 1:size(eleccolors)
          colors{i} = eleccolors(i,:);
          colourIndex(i) = i;
       end
   end
   colourMap = eleccolors; % pass back the input as the output
end

if plot
    
    mriColor = [255 224 179]/256;
    openepColor = [.7 .7 .7];
    
    % atrium with barycentre lines
    figure
    drawMap(userdata, 'type', 'none', 'usrcolormap', openepColor);
    hold on
    for i = 1:10:size(vertices,1)
        line( [C_ep(1,1) vertices(i,1)], [C_ep(1,2) vertices(i,2)], [C_ep(1,3) vertices(i,3)], 'linewidth', 1.5, 'color', 'k');
        vert4points(i,:) = vertices(i,:);
    end
    vert4points(vert4points(:,1)==0,:) = [];
    plot3(vert4points(:,1), vert4points(:,2), vert4points(:,3), '.', 'color', colorBrewer('r'), 'markersize', 15);
    plotsphere(C_ep(1,1), C_ep(1,2), C_ep(1,3),'k',1,16);
    title('OpenEP geometry with subset of barycentre lines')

    % 3D plot of selected electrodes on OpenEP surface
    figure
    hSurf = drawMap(userdata, 'type', 'none', 'usrcolormap', openepColor);
    hold on
    for i = 1:numel(elecsamples)
         coord = X(elecsamples(i),:)/1000;
         plotTag(userdata, 'coord', coord, 'color', colors{i}, 'size', 1.5);
         line( [coord(1) rV(elecsamples(i),1)], [coord(2) rV(elecsamples(i),2)], [coord(3) rV(elecsamples(i),3)], 'linewidth', 1.5, 'color', 'k');
         line( [coord(1) C_ep(1,1)],            [coord(2) C_ep(1,2)],            [coord(3) C_ep(1,3)],            'linewidth', 1.5, 'color', [.3 .3 .3], 'linestyle', '--');
    end
    plotsphere(C_ep(1,1), C_ep(1,2), C_ep(1,3),'k',1,16);
    title('OpenEP geometry with selected electrodes and barycentre')
    set(hSurf, 'facealpha', 0.8);
    
    % UAC plot of selected electrodes - OpenEP surface
    figure
    hSurf = trisurf(getTriangulationFromMeshStruct(openEpUacMesh));
    set(hSurf, 'facecolor', openepColor, 'edgecolor', 'none');
    hold on
    for i = 1:numel(elecsamples)
        plotTag(userdata, 'coord', coordUac(elecsamples(i),:), 'color', colors{i}, 'size', 3/100);
    end
    title('UAC representation of OpenEP geometry with surface-projected electrode positions')
    view(0,90);
    xlabel('UAC-1');ylabel('UAC-2');
    set(gcf, 'color', 'w')
    
    % UAC plot of selected electrodes - OpenEP distances
    figure
    hSurf = trisurf(getTriangulationFromMeshStruct(openEpUacMesh));
    set(hSurf, 'facecolor', openepColor, 'edgecolor', 'none');
    hold on
    for i = 1:numel(elecsamples)
        coord = coordUac(elecsamples(i),:);
        coord(3) = dR(elecsamples(i));
        plotTag(userdata, 'coord', coord, 'color', colors{i}, 'size', 3/100);
        line( [coord(1) coord(1)], [coord(2) coord(2)], [0 dR(elecsamples(i))], 'linewidth', 1.5, 'color', 'k');
        line( [coord(1) coord(1)], [coord(2) coord(2)], [dR(elecsamples(i)) 1], 'linewidth', 1.5, 'color', [.3 .3 .3], 'linestyle', '--');
        
    end
    title('UAC representation of OpenEP geometry with electrodes ')
    xlabel('UAC-1'); ylabel('UAC-2'); zlabel('UAC-3');
    set(gcf, 'color', 'w')
    set(gca, 'zlim', [0,1]);
    
    % UAC plot of selected electrodes - MRI surface
    figure
    hSurf = trisurf(getTriangulationFromMeshStruct(imgUacMesh));
    set(hSurf, 'facecolor', mriColor, 'edgecolor', 'none');
    hold on
    for i = 1:numel(elecsamples)
        plotTag(userdata, 'coord', coordUac(elecsamples(i),:), 'color', colors{i}, 'size', 3/100);
    end
    title('UAC representation of MRI geometry with surface-projected electrode positions')
    view(0,90);
    xlabel('UAC-1');ylabel('UAC-2');
    set(gcf, 'color', 'w')
    
    % UAC plot of selected electrodes - MRI distances
    figure
    hSurf = trisurf(getTriangulationFromMeshStruct(imgUacMesh));
    set(hSurf, 'facecolor', mriColor, 'edgecolor', 'none');
    hold on
    for i = 1:numel(elecsamples)
        coord = coordUac(elecsamples(i),:);
        coord(3) = dR(elecsamples(i));
        plotTag(userdata, 'coord', coord, 'color', colors{i}, 'size', 3/100);
        line( [coord(1) coord(1)], [coord(2) coord(2)], [0 dR(elecsamples(i))], 'linewidth', 1.5, 'color', 'k');
        line( [coord(1) coord(1)], [coord(2) coord(2)], [dR(elecsamples(i)) 1], 'linewidth', 1.5, 'color', [.3 .3 .3], 'linestyle', '--');
        
    end
    title('UAC representation of MRI geometry with electrodes ')
    xlabel('UAC-1'); ylabel('UAC-2'); zlabel('UAC-3');
    set(gcf, 'color', 'w')
    set(gca, 'zlim', [0,1]);
    
    % 3D plot of selected electrodes on MRI surface
    figure
    hSurf = trisurf(getTriangulationFromMeshStruct(img3DMesh));
    hold on
    drawFreeBoundary(getTriangulationFromMeshStruct(img3DMesh,'type','trirep'),'k');
    set(hSurf, 'facecolor', mriColor, 'edgecolor', 'none');
    cameraLight;
    axis off equal vis3d;
    hold on
    for i = 1:numel(elecsamples)
         coord = Y(elecsamples(i),:);
         plotTag(userdata, 'coord', coord, 'color', colors{i}, 'size', 3000);
         line( [coord(1) nV(elecsamples(i),1)], [coord(2) nV(elecsamples(i),2)], [coord(3) nV(elecsamples(i),3)], 'linewidth', 1.5, 'color', 'k');
         line( [coord(1) C_img(1,1)],           [coord(2) C_img(1,2)],           [coord(3) C_img(1,3)],            'linewidth', 1.5, 'color', [.3 .3 .3], 'linestyle', '--');
    end
    plotsphere(C_img(1,1), C_img(1,2), C_img(1,3),'k',1000,16);
    title('OpenEP geometry with selected electrodes and barycentre')
    set(gcf, 'color', 'w');
end

end