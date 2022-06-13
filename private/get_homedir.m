function hd = get_homedir()
%GET_HOMEDIR Returns the user's home directory.

if ispc
    hd = [getenv('HOMEDRIVE') getenv('HOMEPATH')];
else
    hd = getenv('HOME');
end

end