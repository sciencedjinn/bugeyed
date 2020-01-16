function r=loadmov(FOLDERNAME)
%LOADMOV(Foldername) loads a movie from bmp-files
%Movie must be in work\Movies\Foldername

if ~exist([matlabroot, '\work\Movies', FOLDERNAME])
    error(['*** Error: Folder ', FOLDERNAME, ' does not exist']);
end

disp(['* Loading movie ', FOLDERNAME, ' ...']);
i=0;
currentFilename=[matlabroot, '\work\Movies', FOLDERNAME, '\', FOLDERNAME, int2str(i), '.bmp'];
while exist(currentFilename)
    currentFile=imread(currentFilename);                                    %Read current input file
    X(i+1)=im2frame(currentFile);                                             %Concatenate new File to InputData
    i=i+1;                                                                  %Go to next filenumber
    currentFilename=[matlabroot, '\work\Movies', FOLDERNAME, '\', FOLDERNAME, int2str(i), '.bmp'];      %And calculate next filename
end
disp([sprintf('\b'), 'done']);
r=X;    %return output