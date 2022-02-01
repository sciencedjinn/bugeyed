function bugeyed(inGeomName, inVisName, inPara)
    % BUGEYED 1.01 filters a video or single image to simulate arthropod vision.
    % Data of the eye geometry have to be given in geometrical input files in asciiformat in be_input\Geom, delimited by commas.
    % Visual input has to be stored as image sequences in be_input\vis.
    % For further information see default definitons below.
    %
    % Example:
    % bugeyed('Crab', 'Bird_0')
    %
    %    inGeomName - can be a full path (at least two directories deep) or a pre-defined geometry name (must be found in be_input/geom)
    %    inVisName  - can be a full path (at least two directories deep) or a pre-defined stimulus name (must be found in be_input/vis)


    bugeyed_paths;

    %% defaults
    if nargin<1 || isempty(inGeomName), inGeomName = sub_getInputName('geom'); end % ask user
    if nargin<2 || isempty(inVisName),  inVisName  = sub_getInputName('vis'); end % ask user
    if nargin<3 || isempty(inPara)
        inPara = bugeyed_defaultPara;
    end

    %% load eye data
    fprintf('\nCalculating input image/video ''%s'' as seen by eye geometry ''%s''\n\n', inVisName, inGeomName);

    fprintf('* Loading eye geometry...\n');
    geom = sub_loadGeom(inGeomName, inPara); % load optical axes and acceptance angles
    fprintf('\bdone.\n');

    if inPara.stimulusMode
        %% load and filter visual input files
        visFiles = sub_findVisualInputFiles(inVisName, inPara.imFormat);  % all image file names
        visNumber = length(visFiles);                                     % total number of image files
        if visNumber==0, error('No files of %s type found in folder %s', inPara.imFormat, inVisName); end
        if visNumber>1, error('Stimulus mode only supports a single image!'); end
        
        tic
        im = sub_loadVisInput(visFiles, inPara.verbose);
        if inPara.gradientOnly
            im = repmat(mean(im, 1), [size(im, 1) 1]);
        end

        t = 0:1/inPara.frameRate:inPara.stimulus.duration;
        ims = nan(size(im, 1), size(im, 2), 3, length(t)); % pre-allocate
        for i = 1:length(t)
            ims(:, :, :, i) = bugeyed_overlayStimulus(im, inPara.visualField, inPara.stimulus, t(i));
        end

        % Filter
        fprintf('* Spatially filtering input files...\n');
        FilteredVisData = sub_filterSpatially(geom, ims, inPara);
        fprintf('\b done.\n');
        visNumber = length(t);
    else
        %% load and filter visual input files
        visFiles = sub_findVisualInputFiles(inVisName, inPara.imFormat);  % all image file names
        visNumber = length(visFiles);                                     % total number of image files
        if visNumber==0, error('No files of %s type found in folder %s', inPara.imFormat, inVisName); end
    
        a = 1;
        tic
        while a <= visNumber
            b = min([a + inPara.bufferSize - 1, visNumber]);
    
            % Load
            fprintf('* Loading visual input files %d to %d out of %d...\n', a, b, visNumber);
            ims = sub_loadVisInput(visFiles(a:b), inPara.verbose);
            
            if inPara.gradientOnly
                ims = repmat(mean(ims, 1), [size(ims, 1) 1 1]);
            end
    
            fprintf('\b done.\n');
    
            if a == 1
                FilteredVisData = nan(size(geom.OAs, 1), size(ims, 3), length(visFiles)); % pre-allocateonce you know the channel number
            end
    
            % Filter
            fprintf('* Spatially filtering input files %d to %d out of %d...\n', a, b, visNumber);
            FilteredVisData(:, :, a:b) = sub_filterSpatially(geom, ims, inPara);
            fprintf('\b done.\n');
    
            a = b + 1;
        end
    end

    %% plot and save

    % Prepare figure
    figure(4); 
    clf; 
    h = axes; 
    maxfig;

    if size(FilteredVisData, 2) == 1
        colormap gray;
    end

    if ~exist(inPara.saveDir, 'file')
        mkdir(inPara.saveDir);
    end

    [~, shortGeomName] = fileparts(inGeomName);
    [~, shortVisName]  = fileparts(inVisName);
    if visNumber==1
        % Just a single image, write to image file
        saveName = fullfile(inPara.saveDir, [shortGeomName, 'View', shortVisName, inPara.outImageFormat]);
        
        showVoronoi(inGeomName, inPara, FilteredVisData(:, :, 1)./inPara.normalisation, h, [], inPara.forceRecalc);
        if inPara.plotType=='p', axis(h, inPara.visualField); end
        imwrite(getframe(h).cdata, saveName);
    else
        % Write to output video
        saveName = fullfile(inPara.saveDir, [shortGeomName, 'View', shortVisName, inPara.outVideoFormat]);
        v = VideoWriter(saveName, 'MPEG-4');
        v.FrameRate = inPara.frameRate;
        open(v);
        ph = []; % empty handle means it will be created the first time showVoronoi is called
        for i = 1:visNumber
            ph = showVoronoi(inGeomName, inPara, FilteredVisData(:, :, i)./inPara.normalisation, h, ph, inPara.forceRecalc&&i==1);
            title(sprintf('Image #%3d', i));
            drawnow;
            axis(h, 'equal')
            if inPara.plotType=='p', axis(h, inPara.visualField); end
            writeVideo(v, getframe(h)); % saving frames in an array and moving the writeVideo out of the loop only saves a few seconds - not worth it
        end
        close(v);
    end

end % main






%% SUBFUNCTIONS

%% SUB_GETINPUTNAME
function outDefaultName = sub_getInputName(inType)
    % Poll the user with a choice of input files
    Names = struct2cell(dir(fullfile(bugeyed_rootFolder, 'be_input', inType)));
    SelectionList = squeeze(Names(1, 3:size(Names, 2)));
    switch inType
        case 'vis'
            SelectionText = 'Select image sequence';
        case 'geom'
            SelectionText = 'Select geometry';
    end
    [InputNumber, ok] = listdlg('ListString', SelectionList, 'ListSize', [300, 300], 'SelectionMode', 'single', 'Name', SelectionText, 'CancelString', 'Cancel');
    
    % routine for "Cancel..."
    if ok == 0
        error('To use data not included here, please use the ''Import''-function.');
    else
        outDefaultName = SelectionList{InputNumber};
    end               
    close all;
end

%% SUB_FINDVISUALINPUTFILES
function outFilenames = sub_findVisualInputFiles(inVisName, visFormat)
    % Find all image files with an appropriate extension

    if isfolder(inVisName)
        fnMask = fullfile(inVisName, ['*', visFormat]);
        outFilenames = arrayfun(@(x) fullfile(x.folder, x.name), dir(fnMask), 'UniformOutput', false);
    elseif isfile(inVisName)
        outFilenames = {inVisName};
    else
        % old standard behaviour, look in the local stimulus library
        visPath = fullfile(bugeyed_rootFolder, 'be_input', 'vis', inVisName);
        fnMask = fullfile(visPath, ['*', visFormat]);
        outFilenames = arrayfun(@(x) fullfile(x.folder, x.name), dir(fnMask), 'UniformOutput', false);
    end
end                          

%% SUB_LOADVISINPUT
function ims = sub_loadVisInput(fileNames, verbose)                             
    % Subroutine 'Visual Input', loads visual input files and concatenates them
    if verbose, h = waitbar(0, 'Loading... (0%%)', 'Name', 'Progress'); end
    for fileNumber = 1:length(fileNames)
        im = double(imread(fileNames{fileNumber})); % read current input file
        im = fliplr(permute(im, [2 1 3]));                           % get indices in correct order (so azimuth is still first index)
        
        if fileNumber == 1
            ims = zeros([size(im, 1), size(im, 2), size(im, 3), length(fileNames)]); % pre-allocate
        end
        ims(:, :, :, fileNumber) = im;                                 % Concatenate new File to InputData
        if verbose
            perc = fileNumber/length(fileNames);
            waitbar(perc, h, sprintf('Loading... (%d%%, %ds remaining)', round(100*perc), round(toc/perc-toc)));
        end
    end
    if verbose, close(h); end
end

%% SUB_LOADGEOM
function outGeom = sub_loadGeom(inGeomName, inPara)
    % Subroutine 'Eye Geometry Input', loads geometrical files, which have to be in a very specific format

    pointsFileName = bugeyed_fileName(inGeomName, 'points');
    anglesFileName = bugeyed_fileName(inGeomName, 'angles');
    diamsFileName  = bugeyed_fileName(inGeomName, 'diams');
    ioasFileName   = bugeyed_fileName(inGeomName, 'ioas');

    outGeom.OAs    = dlmread(pointsFileName); % load optical axes
    
    % Load lens diameter or inter-ommatidial angle data
    if exist(anglesFileName, 'file')
        % if acceptanceAngle data exists, load it and use it preferentially
        % this must be an estimated half-width (FWHM) for each receptor/facet
        outGeom.rho   = dlmread(anglesFileName);   % acceptance angle rho, in degrees (Airy disc FWHM)
    elseif exist(diamsFileName, 'file')
        % lens diameter data has the next highest priority
        outGeom.diams = dlmread(diamsFileName);
        outGeom.rho   = rad2deg(1.02 * inPara.waveLength/1000 ./ outGeom.diams);   % acceptance angle rho, in degrees (Airy disc FWHM)
    elseif exist(ioasFileName, 'file')
        % otherwise, try to load inter-ommatidial angle data
        outGeom.IOAs = dlmread(ioasFileName);
        outGeom.rho   = inPara.phi2rho .* outGeom.IOAs;                         % acceptance angle rho, in degrees (Smolka2004,p.21)
    else
        error('No acceptance angle, lens diameter or IOA data found.');
    end
    % Calculate Gaussian st.d.
    outGeom.sigma = 1/(2*sqrt(2*log(2))) .* outGeom.rho; % Gaussian standard deviation sigma, in degrees
end  
 
%% SUB_FILTERSPATIALLY
function outFilteredInput = sub_filterSpatially(inGeom, inIms, inPara)
    % Subroutine sub_filterSpatially computes a quantum-catch/pixel-catch for each  

    % Parameters
    cNum        = size(inIms, 3);           % number of channels in input image
    fNum        = size(inIms, 4);           % number of timesteps/frames
    azi         = linspace(inPara.visualField(1), inPara.visualField(2), size(inIms, 1)+1); % vector of azimuthal positions in input image pixel edges
    azi         = azi(1:end-1)+diff(azi)/2;                                                 % vector of azimuthal positions in input image pixel centres
    ele         = linspace(inPara.visualField(3), inPara.visualField(4), size(inIms, 2)+1); % vector of elevation positions in input image pixel edges
    ele         = ele(1:end-1)+diff(ele)/2;                                                 % vector of azimuthal positions in input image pixel centres
    dpp_a       = median(abs(diff(azi)));   % degrees per pixel along azimuth of original images
    dpp_e       = median(abs(diff(ele)));   % degrees per pixel along elevation of original images
    Area        = (dpp_e*ones(length(azi), 1)) * (dpp_a*cosd(ele)); % The area (in degrees^2) covered by each pixel in the image
    
    % determine which ommatidia need to be calculated
    sel = ~(inGeom.OAs(:, 1)<inPara.ommatidialField(1) |...
            inGeom.OAs(:, 1)>inPara.ommatidialField(2) |...
            inGeom.OAs(:, 2)<inPara.ommatidialField(3) |...
            inGeom.OAs(:, 2)>inPara.ommatidialField(4)); 
    ommInds = find(sel);
    nRelOmms = length(ommInds);
    
    % pre-allocate/initialise loop variable
    outFilteredInput = zeros(size(inGeom.OAs, 1), cNum, fNum);
    if inPara.verbose
        perc = 0;
        h = waitbar(0, 'Spatially filtering... (0%%)', 'Name', 'Progress');
    end

    try
        for i = 1:nRelOmms        % for each relevant ommatidium
            ommInd = ommInds(i);            % index in inOAs/inIOAs for this ommatidium
            cAz    = inGeom.OAs(ommInd, 1); % optical axis azimuth, in degrees
            cEl    = inGeom.OAs(ommInd, 2); % optical axis elevation, in degrees
            corr   = 1 / cosd(cEl);         % correction factor for elevation in equirectangular projection
            [rows, cols, gk] = nested_calcGaussian;

            for f = 1:fNum % for each frame
                for c = 1:cNum % for each channel
                    V = inIms(cols, rows, c, f); % (Indexing this earlier is slower)
                    outFilteredInput(ommInd, c, f) = sum(sum(V.*gk));  % Filters every timestep with a Gauss filter 
                end
            end
        
            if inPara.verbose
                newPerc = round(100*i/nRelOmms);
                if newPerc>perc
                    perc = newPerc;
                    waitbar(perc/100, h, sprintf('Spatially filtering... (%d%%, %ds remaining)', perc, round(100*toc/perc-toc)));
                end
            end
        end
        if inPara.verbose, close(h); end
    catch me
        if inPara.verbose, close(h); end
        rethrow(me);
    end

    % subfunctions
    function [rows, cols, gk] = nested_calcGaussian
                % cut a 6 * sig square from the image
                minAz = cAz - inPara.filterCutOff * corr * inGeom.sigma(ommInd);
                maxAz = cAz + inPara.filterCutOff * corr * inGeom.sigma(ommInd);
                minEl = cEl - inPara.filterCutOff * inGeom.sigma(ommInd);
                maxEl = cEl + inPara.filterCutOff * inGeom.sigma(ommInd);
                if maxAz<min(azi) || minAz>max(azi) || maxEl<min(ele) || minEl>max(ele)
                    % square lies completely outside the target area, set this receptor to 0
                    rows = 1;
                    cols = 1;
                    gk = 0;
                    return;
                end
    
                lowrow      = sub_findfirstbelow(ele, max([minEl min(ele)]));  % first row below ele - 3 * sigma
                highrow     = sub_findfirstabove(ele, min([maxEl max(ele)]));  % first row above ele + 3 * sigma
                lowcol      = sub_findfirstbelow(azi, max([minAz min(azi)]));  % first col below azi - 3 * sigma * corr
                highcol     = sub_findfirstabove(azi, min([maxAz max(azi)]));  % first col above azi + 3 * sigma * corr
                rows        = min([highrow lowrow]):max([highrow lowrow]);  % numbers of all rows between highrow and lowrow
                cols        = min([highcol lowcol]):max([highcol lowcol]);  % numbers of all rows between highrow and lowrow
                kernelEles  = ele(rows);
                kernelAzis  = azi(cols);
        
                AreaX       = Area(cols, rows);
                [KA, KE]    = ndgrid(kernelAzis, kernelEles);
                gk          = exp( -spdist(cAz, cEl, KA, KE).^2 / (2*inGeom.sigma(ommInd).^2) ); % Gaussian kernel, up to a distance of 3 * sig (corrected for ele)
                gk          = gk .* AreaX;                                                  % weight kernel by relative pixel areas
                gk          = gk / sum(gk(:));                                              % normalise
    end

    function pos = sub_findfirstbelow(mat, el)
        % finds the elements in mat nearest to but smaller than el
        mat(mat>el) = Inf;
        [~, pos] = min(abs(mat-el));
    end

    function pos = sub_findfirstabove(mat, el)
        % finds the elements in mat nearest to but larger than el
        mat(mat<el) = Inf;
        [~, pos] = min(abs(mat-el));
    end
end
        
