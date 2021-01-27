function [ faceData ] = trVertToFaceData( tr, vertexData )
% TRVERTTOFACEDATA converts vertex to face data for a TriRep object.
% Usage:
%   faceData = trVertToFaceData(tr, vertexData)
% Where:
%   tr - is a TriRep object with the following dimensions:
%           tr.Triangulation    q * 3
%           tr.X                p * 1
%   vertexData - is an p*1 array of scalar data (by vertex)
%   faceData - is a q*1 array of scalar data (by face)
%
% TRVERTTOFACEDATA detailed description.
%
% Author: Steven Williams (2012)
% Modifications - 
%
% Info on Code Testing:
						% ---------------------
                        % test code
                        % ---------------------
%
% Use warnings with the following format
%    disp('TRVERTTOFACEDATA: warning.')
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

%faceData = zeros(length(tr.Triangulation), 1);
% for i = 1:length(tr.Triangulation);
%     v1 = tr.Triangulation(i, 1);
%     v2 = tr.Triangulation(i, 2);
%     v3 = tr.Triangulation(i, 3);
%     
%     faceData(i) = mean([vertexData(v1) vertexData(v2) vertexData(v3)]);
% end

faceData = mean(vertexData(tr.Triangulation),2);