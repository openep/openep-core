function listing = nameFolds(pathname)
% NAMEFOLDS lists all subfolders within pathname.
% Usage:
%   listing = nameFolds(pathname)
% Where:
%   a is the input
%   b is the output
%
% NAMEFOLDS lists all subfolders within a directory by using the matlab
% command dir and removing anything that isn't a directory.
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

d = dir(pathname);
isub = [d(:).isdir]; % returns logical vector
listing = {d(isub).name}';
listing(ismember(listing,{'.','..'})) = [];

end