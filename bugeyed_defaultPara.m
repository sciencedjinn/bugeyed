function bePara = bugeyed_defaultPara()
% BUGEYED_DEFAULTPARA returns to default parameter object for bugeyed plotting
% Any of the parameters can be overwritten before sending them to bugeyed,
% BUT DO NOT OVERWRITE ANY PARAMETERS IN THIS FUNCTION!!!!
%
% Example:
%    Para = bugeyed_defaultPara;
%    Para.RhabdomsToExclude = [1011 1361];
%    Para.RhabdomsToReverse  = [1161 1132 598 1305]; % Flip the direction of these rhabdoms
%    fn = 'Av1h1_half_data';
%    wormResolution (fn, Para);

%% PARA object
bePara.verbose = true;                      % toggle the display of waitbars during loading
bePara.forceRecalc = false;

%% visual input parameters
bePara.visualField = [-180 180 -90 90];     % Camera visual field: [minAzi, maxAzi, minEle, maxEle] (degrees)
bePara.ommatidialField = [-180 180 -90 90]; % Ommatidia with optical axes inside this field are calculated: [minAzi, maxAzi, minEle, maxEle] (degrees)
bePara.imFormat = '.jpg';                   % file format of visual input files
bePara.bufferSize = 100;                    % for image sequence input: number of images to buffer before filtering. Larger number excutes faster, but requires more RAM.
bePara.normalisation = 2^16-1;              % normalise the output video to this maximum (2^8-1 for normal jpgs, 2^16-1 for 16-bit tifs)
                                            % use a higher or lower value to decrease/increase output brightness at the risk of saturation (currently only works for rgb inputs)
bePara.gradientOnly = false;                % Set this to true to remove from the input images any horizontal variation (e.g. for ELF mean images)

%% Stimulus 
% If the input is a single image, stimulusMode can be used to overlay a moving stimulus on this image and create an output video of that
bePara.stimulusMode = false;
bePara.stimulus.type = 'LoomingDot';        % Currently supported: LoomingDot and MovingDot
bePara.stimulus.initialPosition = [-30 0];  % initial stimulus azimuth/elevation position, in degrees
bePara.stimulus.initialDiameter = 10;       % initial stimulus diameter (degrees)
bePara.stimulus.movementSpeed = [20 0];     % for MovingDot ONLY; stimulus movement speed along azimuth and elevation (degrees/s); 
                                            % Warning: If both are set (diagonal movement), the result will not be geometrically meaningful!
bePara.stimulus.timeToContact = 4;          % for LoomingDot ONLY; time to collision, in seconds
bePara.stimulus.duration = 4;               % duration of the stimulus movement, in seconds
bePara.stimulus.opacity = 1;                % opacity of the stimulus, in log10; e.g. opacity 2 means the stimulus location will be dimmed by a factor of 100

%% Acceptance angles
% Receptor acceptance angles can be provided directedly (in a *Angles.txt file) or calculated from inter-ommetidial angles (*IOA.txt) or facet diameters (*Diam.txt)
bePara.phi2rho = 1;                         % assumed ratio between inter-ommatidial angle (IOA = phi) and receptor acceptance angle (AA = rho).
                                            % Influences acceptance angles ONLY if calculated from inter-ommatidial angles.
bePara.waveLength = 450;                    % assumed wavelength of processed light. Influences acceptance angles ONLY if calculated from facet diameter.
bePara.filterCutOff = 3;                    % acceptance angle filter kernels will be cut off after this many standard deviations (default 3 to include 99.7% of data)

%% Plotting
bePara.maxVoronoiDist = 180;
bePara.plotView = [60, 20];                 % The azimuth and elevation from which the 3d plot is viewed
bePara.plotType = 's';                      % How to plot the output images: 's' for sphere (3d plot), 'p' for plane (2d plot);
bePara.voronoiEdgeLimit = 60;               % Limit for edge length to plot in the voronoi diagram; this is necessary to limit extreme-elevation 
                                            % facets from being plotted in 2d mode (they will not come out right no matter what)

%% visual output parameters
bePara.frameRate = 25;                      % output video frame Rate (Hz)
bePara.saveDir = fullfile(bugeyed_rootFolder, 'output'); % output files will be saved in this folder
bePara.outImageFormat = '.jpg';             % file format of output images
bePara.outVideoFormat = '.mp4';             % file format of output videos


