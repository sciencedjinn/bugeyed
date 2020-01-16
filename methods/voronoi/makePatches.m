function makePatches(inGeomName)

vFolder = fullfile(rootFolder, 'be_input', 'voronoi', inGeomName);
vFileIn = fullfile(vFolder, [inGeomName, 'Voronoi.mat']);
vFileOut = fullfile(vFolder, [inGeomName, 'VoronoiPatches.mat']);

%% load Voronoi and convert to a quick-to-plot format
load(vFileIn, 'outVOR');

polSize = cellfun(@(x) size(x, 1), outVOR); % number of vertices for each polygon
for np = min(polSize):max(polSize)
    % find all polygons with np vertices
    theseVerts = find(polSize==np);

    Fx{np} = cell2mat(cellfun(@(x) x(:, 1), outVOR(theseVerts), 'UniformOutput', false)); %#ok<*AGROW>
    Fy{np} = cell2mat(cellfun(@(x) x(:, 2), outVOR(theseVerts), 'UniformOutput', false));
    Id{np} = theseVerts;
end

save(vFileOut, 'Fx', 'Fy', 'Id')