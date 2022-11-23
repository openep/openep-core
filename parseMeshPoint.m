function position = parseMeshPoint(unityFile, pointName)
% PARSEMESHPOINT Extracts the locations of mesh points from the Unity.txt
% file of a Kodex export dataset
%
% Usage:
%   position = parseMeshPoint(unityFile, pointName)
% Where:
%   unityFile  - absolute path to the Unity.txt file
%   pointName  - the required point name
%   position  - the point co-ordinates
%
% PARSEMESHPOINT does not accept any parameter-value pairs
%
% PARSEMESHPOINT is required when mesh points have been added manually
% during a Kodex case, for example to define a pacing site on the mesh.
% These points are not saved as landmark points but their locations are
% encoded within the Unity.txt. A two-step process is needed to extract
% their co-ordinates. First, the point_key has to be identified. Secondly,
% the position of the point labelled point_key has to be identified. Note
% that the y and z co-ordinates of the position need to be swapped to
% correctly locate the point in the Kodex mesh space.
%
% Author: Steven Williams (2022)
% Modifications -
%
% Info on Code Testing:
% ---------------------------------------------------------------
% unityFile = '/Users/Steven/Desktop/2021-08-03_12-59-54_PL_RA_LAT_WV_ABL/task_manager_logs/20210803_0906_26/Unity.txt';
% pointName = 'CS56_P500';
% position = parseMeshPoint(unityFile, pointName)
% ---------------------------------------------------------------
%
% ---------------------------------------------------------------
% code
% ---------------------------------------------------------------

fid = fopen(unityFile);
lineNum = 0;
tLine = {};
while 1
    lineNum = lineNum+1;
    tLine{lineNum} = fgetl(fid);
    if ~ischar(tLine{lineNum})
        break
    end
end
fclose(fid);
tLine = tLine';

ispresent = cellfun(@(s) ~isempty(strfind(s, pointName)), tLine); %#ok<*STREMP>
iLines = find(ispresent);

if isempty(iLines)
    % The requested point name was not found in the Unity.txt file so
    % return an error
    error(['OPENEP/PARSEMESHPOINT.M: A point named ' pointName ' was not found in the Unity.txt file']);
else
    % The requested point name was found in the Unity.txt file so look for
    % its co-ordinates

    % Begin by finding the point_key for every line where the point name
    % was found
    for i = 1:numel(iLines)

        thisLine = tLine{iLines(i)};
        iStart = strfind(thisLine, 'point_key') - 1;
        iSquareClose = strfind(thisLine, ']');
        iEnd = iSquareClose(1) - 2;
        requiredString = thisLine(iStart:iEnd);
        stringComponents = strsplit(requiredString, ',');

        % by definition, point_key will appear in the first cell of
        % stringComponents, use regexp to pull out the point_key
        components = regexp(stringComponents{1}, '"(.*?)"', 'tokens');
        point_key{i} = components{2}{1};

    end
    point_key = unique(point_key);
end

if numel(point_key) > 1
    % We have not handled the situation where there are identical points
    % labelled with the same name yet
    error('OPENEP/PARSEMESHPOINT.M: Multiple points found with the same name is not handled yet');
else
    % Find all lines referencing the required point
    ispresent = cellfun(@(s) ~isempty(strfind(s, point_key{:})), tLine);
    iLines = find(ispresent);

    % Find all lines with position data
    ispresent = cellfun(@(s) ~isempty(strfind(s, 'position')), tLine);
    iPosition = find(ispresent);

    % Find the line for this point which has position data
    iLine = intersect(iLines, iPosition);
    thisLine = tLine{iLine};

    % Parse the position
    iStart = strfind(thisLine, 'position') - 1;
    iSquareClose = strfind(thisLine, ']');
    iEnd = iSquareClose(2) - 1; % use index (2) this time since we want to include the first square bracket occuring after the position co-ordinates
    requiredString = thisLine(iStart:iEnd);

    coords = regexp(requiredString, '\[(.*?)\]','tokens');
    coords = strsplit(coords{1}{1},',');

    x = str2double(coords{1});
    y = str2double(coords{2});
    z = str2double(coords{3});

    % switch the y and z co-ordinates
    position = [x z y];

end

end