function dropboxDir = localDropboxPath()
    % Get dropboxDir & intermediate data files location
    computerInfo = GetComputerInfo();
    switch (computerInfo.localHostName)
        case 'Ithaka'
            dropboxDir = '/Volumes/SSDdisk/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        case 'Crete'
            dropboxDir = '/Volumes/Dropbox/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        case 'Santorini'
            dropboxDir = '/Users/nicolas/Aguirre-Brainard Lab Dropbox/Nicolas Cottaris';

        otherwise
            if (contains(computerInfo.networkName, 'leviathan'))
                dropboxDir = '/media/dropbox_disk/Aguirre-Brainard Lab Dropbox/isetbio isetbio';
            else
                error('Could not establish dropbox location')
            end
    end
end