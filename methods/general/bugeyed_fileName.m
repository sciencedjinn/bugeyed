function [fileName, filePath, shortGeomName] = bugeyed_fileName(inGeomName, fileType)
% returns the appropriate filenames and paths for geometry files
% submit a relative or full path (at least two directories deep) to follow that; otherwise, searches in the built-in geom database

[filePath, shortGeomName] = fileparts(inGeomName);

switch lower(fileType)
    case {'points', 'angles', 'diams', 'ioas'}
        if isempty(filePath)
            % old standard behaviour, look in the local geometry library
            filePath = fullfile(bugeyed_rootFolder, 'be_input', 'geom', shortGeomName);
        end
        switch lower(fileType)
            case 'points'
                fileName = fullfile(filePath, [shortGeomName, 'Points.txt']);
            case 'angles'
                fileName = fullfile(filePath, [shortGeomName, 'Angles.txt']);
            case 'diams'
                fileName = fullfile(filePath, [shortGeomName, 'Diams.txt']);
            case 'ioas'
                fileName = fullfile(filePath, [shortGeomName, 'IOAs.txt']);
        end

    case {'v', 'voronoi', 'vor', 'p', 'patch', 'patches'}
        if isempty(filePath)
            % old standard behaviour, look in the local geometry library
            filePath = fullfile(bugeyed_rootFolder, 'be_input', 'voronoi', shortGeomName);
        end
        switch lower(fileType)
            case {'v', 'voronoi', 'vor'}
                fileName = fullfile(filePath, [shortGeomName, 'Voronoi.mat']);
            case {'p', 'patch', 'patches'}
                fileName = fullfile(filePath, [shortGeomName, 'VoronoiPatches.mat']);
        end
    otherwise
        error('Unknown file type: %s', fileType);
end
