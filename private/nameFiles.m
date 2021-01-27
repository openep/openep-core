function listing = nameFiles(pathname, varargin)
% NAMEFFILES lists all files within pathname.
% Usage:
%   listing = nameFolds(pathname)
% Where:
%   listing - is a cell array of filenames
%   pathname - is the directory in question
%
% NAMEFILES lists all files within a directory by using the matlab
% command dir and removing anything that isn't a directory.
%
% Parameter/Value pairs:
%   'showhiddenfiles' {true}|false
%   'extension' string
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

%defaults
showHiddenFiles = true;

%parse inputs
p = inputParser;
p.addRequired('pathname', @ischar);
p.addParamValue('showhiddenfiles', true, @islogical);
p.addParamValue('extension', [], @ischar);
p.parse(pathname, varargin{:});
inputs = p.Results;
pathname = inputs.pathname;
extension = inputs.extension;

d = dir(pathname);
isub = [d(:).isdir]; % returns logical vector
listing = {d(~isub).name}';
listing(ismember(listing,{'.','..'})) = [];

if ~showHiddenFiles
    listing(strstartcmpi('.', listing)) = [];
end

if ~isempty(extension)
    listing(~strstartcmpi(fliplr(extension), cellfun(@fliplr, listing, 'uniformoutput', false))) = [];
end
