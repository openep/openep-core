function [Y, X] = translateOpenEP2ImagingPositions(userdata, mOpenEP3D, mOpenEpUac, mImg3D, mImgUac, varargin)
% TRANSLATEOPENEP2IMAGINGPOSITIONS Translates electrode positions from OpenEP
% space to Imaging space
%
% Usage:
%   [Y, surfY] = translateOpenEP2ImagingPositions(userdata, openEpUac, img3D, imgUac)
% Where:
%   userdata  - an OpenEP dataset
%   mOpenEP3D - on OpenEP geometry, presented as a CarpMesh object
%   openEPUac - an OpenEP geometry converted into UACs (universal atrial
%               co-ordinates), presented as a CarpMesh object
%   img3D     - an imaging geometry, presented as a CarpMesh object
%   imgUac    - an imaging geometry converted into UACS, presented as a
%               CarpMesh object
%   Y         - 3D Cartesian co-ordinates of the electrodes in imaging space
%   X         - 3D Cartesian co-ordinates of the electrodes in EP space
%
% TRANSLATEOPENEP2IMAGINGPOSITIONS accepts the following parameter-value pairs
%   'plot'      {false} | true
%
% TRANSLATEOPENEP2IMAGINGPOSITIONS Uses the Universal Atrial Co-ordinates
% to convert electrode positions between OpenEP and Imaging co-ordinate
% systems. For example; this function can be used to create surface and 3D
% co-ordinates representing all of the recording electrode positions during
% a clinical case, relative to an imaging dataset such as MRI or CT
%
% Author: Steven Williams (2021)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
%  load('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Carto/Study_1_07_02_2017_19-42-07_1-Map.mat');
%  mOpenEp3D = io_readCARPMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Carto/Labelled');
%  mOpenEpUac = io_readCARPMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Carto/Labelled_Coords_2D_Rescaling_v3_C');
%  mImg3D = io_readCARPMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Model_16/Labelled');
%  mImgUac = io_readCARPMesh('/media/stw11/Data/StevenModellingData/Models/Carto_153_NR/Model_16/Labelled_Coords_2D_Rescaling_v3_C');
%  [X, surfX] = translateOpenEP2ImagingPositions(userdata, mOpenEp3D, mOpenEpUac, mImg3D, mImgUac)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

nStandardArgs = 5;
plot = false;
if nargin > nStandardArgs
    for i = 1:2:nargin-nStandardArgs
        switch varargin{i}
            case 'plot'
                plot = varargin{i+1};
        end
    end
end

% 1. Load the X Cartesian co-ordinates of the electrode positions from
% OpenEP data. Identify the co-ordinates that are internal to the imaging mesh
X_orig = getElectrogramX(userdata);
tr = getTriangulationFromMeshStruct(mOpenEP3D, 'region', 11, 'scale', 'um');
iPoint = getMappingPointsWithinMesh(userdata, 'mask', tr);
X = X_orig(iPoint,:); % X is now all the points of interest but DOES NOT index into userdata.electric

% 2. Find the barycentre lines; these are the set of lines linking every
%    node on the OpenEP-3D mesh to the barycentre of the OpenEP-3D mesh
baryCentre = getCentreOfMass(userdata);
tr = getOpenEPSubsetFromImgMap(userdata, mOpenEP3D, 11);
vertices = tr.X;

%TODO: make the above into a subfunction; returning a node/element
%representation of the barycentre lines

% 3. For every point of interest X; find the barycentre line that it is
%    closest to; and the percentage distance along this line from the
%    OpenEP surface to the barycentre that the point lies
if plot
    hSurf = drawMap(userdata, 'type', 'none');
    set(hSurf, 'facealpha', 0.5);
end

for i = 1:length(X)
    distance = NaN(length(vertices),1);
    for j = 1:length(vertices)
        distance(j) = point_to_line_segment_distance(X(i,:), vertices(j,:), baryCentre);
    end
    [~, iVertex4Electrode(i)] = nanmin(distance);
    
    if plot
        h(1) = plotTag(userdata, 'coord', X(i,:), 'color', 'g'); % plot the electrode recording position
        h(2) = plotTag(userdata, 'coord', baryCentre, 'size', 2);% plot the barycentre
        h(3) = plotTag(userdata, 'coord', vertices(iVertex4Electrode(i),:), 'size', 2); % plot the identified surface point
        % draw a line
        h(4) = line( [vertices(iVertex4Electrode(i),1), baryCentre(1)], [vertices(iVertex4Electrode(i),2), baryCentre(2)], [vertices(iVertex4Electrode(i),3), baryCentre(3)] ...
            , 'color', 'k' ...
            , 'linewidth', 3);
        pause; delete(h);
    end
    
    lineSegmentLength(i) = lineLength([vertices(iVertex4Electrode(i),:); baryCentre]);
    surf2ElectrodeDistance(i) = lineLength(    [vertices(iVertex4Electrode(i),:); X(i,:)]    );
    
    relevantVertex(i,1:3) = vertices(iVertex4Electrode(i),:);
    distanceRatio(i) = surf2ElectrodeDistance(i) / lineSegmentLength(i);
end

% 4. Find the node in the index of the OpenEP-3D point that is attached to
%    the barycentre line of intereste
iVertex_mOpenEP3D = findclosestvertex(mOpenEP3D.Pts/1000, relevantVertex);

% 5. Find the 2D co-ordinates of this node within the OpenEP-UAC mesh
coord_mOpenEPUac = mOpenEpUac.Pts(iVertex_mOpenEP3D,:);

% 6. Find the index of the closest node to this location in the MRI-UAC
%    mesh; this defines the barycentre line of interest
iVertex_mImgUac = findclosestvertex(mImgUac.Pts, coord_mOpenEPUac);

% 7. Calculate the barycentre of the MRI-3D mesh
tempUserdata.surface.triRep = getTriangulationFromMeshStruct(mImg3D, 'type', 'trirep');
C = getCentreOfMass(tempUserdata);

% 8. calcualte vectors from each node toward the barycentre
baryVectors = repmat(C, size(mImg3D.Pts(iVertex_mImgUac))) - mImg3D.Pts(iVertex_mImgUac,:);

% 9. calculate the length of these vectors
baryVectorLengths = vecnorm(baryVectors,2,2);

% 10. create unit vectors for use in calculation
unitVectors = baryVectors./baryVectorLengths;

% 11 Calculate the line between the node of interest and the barycentre of
%    the MRI-3D mesh.
scaledVectors = unitVectors .* repmat(distanceRatio', 1, 3) .* baryVectorLengths;

Y = mImg3D.Pts(iVertex_mImgUac,:) + scaledVectors;

X = X * 1000;
    
    
end