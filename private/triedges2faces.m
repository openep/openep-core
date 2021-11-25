function iFace = triedges2faces(t, e)
% Where
% iFace - indexes into t.ConnectivityList
% t - a triangulation object
% e - an nx2 array of edges, indexing into vertices in t.Points

% make sure edges are all in ascending order
e = sort(e, 2);

% explode the Connectivity List into individual edges
eAll = t.ConnectivityList;
eAll_A = sort(eAll(:, [1 2]), 2);
eAll_B = sort(eAll(:, [2 3]), 2);
eAll_C = sort(eAll(:, [3 1]), 2);

% check if edges in e appear in eAll_A, _B or _C
[~, iA] = ismember(eAll_A, e, 'rows');
[~, iB] = ismember(eAll_B, e, 'rows');
[~, iC] = ismember(eAll_C, e, 'rows');

iFace = iA | iB | iC;

end