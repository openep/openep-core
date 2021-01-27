function tNew = addtrirep(varargin)
%ADDTRIREP creates a new TriRep that contains tA and tB.
% Usage:
%   tNew = addtrirep(tA, tB, tC ....)
%   tNew = addtrirep({tA, tB, tC ...})
% Where:
%   tA,tB etc are TriRep
%
% ADDTRIREP creates a new TriRep that contains tA and tB.
%
% Author: Nick Linton (2011)
% Modifications - 

    if nargin == 1
        tCell = varargin{1};
    else
        tCell = varargin;
    end
    
    nV = 0;
    nT = 0;
    for i = 1:length(tCell)
        if ~isa(tCell{i},'TriRep')
            error('ADDTRIREP: the input must be TriRep objects or a cell array of TriRep objects.')
        end
        nV = nV + size(tCell{i}.X,1);
        nT = nT + size(tCell{i}.Triangulation,1);
    end
    
    newT = zeros(nT,3);
    newX = zeros(nV,3);
    
    currentV = 1;
    currentT = 1;
    for i = 1:length(tCell)
        if ~isa(tCell{i},'TriRep')
            error('ADDTRIREP: the input must be TriRep objects or a cell array of TriRep objects.')
        end
        nV = size(tCell{i}.X,1);
        nT = size(tCell{i}.Triangulation,1);
        newX(currentV:(currentV+nV-1),:) = tCell{i}.X;
        newT(currentT:(currentT+nT-1),:) = tCell{i}.Triangulation + currentV-1;
        currentV = currentV+nV;
        currentT = currentT+nT;
    end
    
    tNew = TriRep(newT, newX);

end



