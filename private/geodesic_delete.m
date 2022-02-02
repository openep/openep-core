% delete mesh, algorithm, or everything at once
% Danil Kirsanov, 09/2007 

function object = geodesic_delete(object)

global geodesic_library;

if ismac
    geodesic_library = 'geodesic_matlab_api_macos';
elseif isunix
    geodesic_library = 'geodesic_matlab_api';
else
    disp('Platform not supported')
end

if ~libisloaded(geodesic_library)       %everything is already cleared
    object = [];
    return;
end

if nargin == 0          % the simplest way to delete everything is to unload the library
   unloadlibrary(geodesic_library);
else
    if libisloaded(geodesic_library)
        if strcmp(object.object_type, 'mesh')       
            calllib(geodesic_library, 'delete_mesh', object.id);      % delete mesh and corresponding algorithms
        else                                        
            calllib(geodesic_library, 'delete_algorithm', object.id); % delete a single algorithm
        end
    end
end
object = [];
