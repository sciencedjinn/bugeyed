function savemov(FOLDERNAME, mov)
%SAVEMOV(Foldername) saves a movie as bmp-files in work\Movies\Foldername

if exist([matlabroot, '\work\Movies\', FOLDERNAME])
    button = questdlg(['MovieSave ', FOLDERNAME, ' already exists. Overwrite?'], 'Movie already exists','Yes','No','No');
    if strcmp(button, 'No')
        disp('** Warning: Movie was NOT saved. Process aborted by user');
        return;
    end
end

disp(['* Saving movie ', FOLDERNAME, ' ...']);
mkdir([matlabroot, '\work\Movies'], FOLDERNAME);
for i=0:size(mov, 2)-1
    [X,Map] = frame2im(mov(i+1));
    currentFilename=[matlabroot, '\work\Movies', FOLDERNAME, '\', FOLDERNAME, int2str(i), '.bmp'];
    imwrite(uint8(X), currentFilename);
end
disp([sprintf('\b'), 'done']);