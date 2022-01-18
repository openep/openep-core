function initialiseGeodesic()
% INITIALISEGEODESIC Initialise the geodesic library and algorithm

%
% Usage:
%   initialiseGeodesic()
%
% Author: Paul Smith (2022)
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

% INITIALISEGEODESIC Initialise the geodesic library and algorithm

global geodesic_library

if ismac
    geodesic_library = 'geodesic_matlab_api_macos';
elseif isunix
    geodesic_library = 'geodesic_matlab_api';
else
    disp('Platform not supported')
end

if ~libisloaded(geodesic_library)
    hfile = 'geodesic_matlab_api.h';
    loadlibrary([geodesic_library '.so'], hfile);
end
