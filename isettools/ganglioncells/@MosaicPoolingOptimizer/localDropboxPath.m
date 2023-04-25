function dropboxDir = localDropboxPath()
    % Get dropboxDir & intermediate data files location
    p = getpref('isetbio');

    if (isfield(p, 'rgcResources'))
        if (isstruct(p.rgcResources))
            if (isfield(p.rgcResources, 'method'))
                switch (p.rgcResources.method)
                    case 'localFile'
                        if (isfield(p.rgcResources, 'URLpath'))
                            dropboxDir = p.rgcResources.URLpath;
                        else
                            error('ISETBio preferences contains an ''rgcResources'' struct, but this struct does not contain a ''URLpath'' field. Check your preferences\n');
                        end
                    otherwise
                        error('Unknown rgcResources.method value: ''%s''. Check your ISETBio preferences.', p.rgcResources.method);
                end
            else
                error('ISETBio preferences contains an ''rgcResources'' struct, but this struct does not contain a ''method'' field. Check your preferences\n');
            end

        else
            error('ISETBio preferences contains an ''rgcResources'' field, but this field is not a struct. Check your preferences\n');
        end

    else
        error('ISETBio preferences does not contain an ''rgcResources'' field. Check your preferences\n');
    end

end