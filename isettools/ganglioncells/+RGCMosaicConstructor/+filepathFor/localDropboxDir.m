function dropboxDirPath = localDropboxDir()

    error('in localDropboxDir() function. Use p = getpref(''isetbio''); resourceDir = p.rgcResources;.')

    dbJsonConfigFile = '~/.dropbox/info.json';
    fid = fopen(dbJsonConfigFile);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);
    dropboxDirPath = val.business.path;
end