function il = illuminantModernize(illuminant)
% Convert old format illuminant structures to the modern format
%
%   illuminant = illuminantModernize(illuminant)
%
% If the structure is already modern, move on.
%
% If it is the old format, it will have the data field wrong.  So, we
% convert here.
% 
% The old format mave these fields:
%
%           data: [31x1 double]  (energy units)
%     wavelength: [1x31 double]
%
% If it also has these string fields, they will be copied too.
%
%        comment: 
%        name:    
%
% Over time, we should convert the data files.  This will be a JEF thing.
%
% (c) Imageval Consulting, LLC 2012

if isfield(illuminant,'type'), 
    % Has a type so, it is modern. Yes, we could test more
    il = illuminant;
    return; 
else
    if ~isfield(illuminant,'data') || ~isfield(illuminant,'wavelength')
        error('No illuminant data or wavelength');
    else
        il = illuminantCreate;
        il = illuminantSet(il,'wavelength',illuminant.wavelength);
        il = illuminantSet(il,'energy',illuminant.data);
        
        if isfield(illuminant,'name')
            il = illuminantSet(il,'name',illuminant.name);
        else
            % This is a fantastic way to name an unknown illuminant.  Find
            % its correlated color temperature and name it that.
            w = illuminantGet(il,'wave');
            spd = illuminantGet(il,'energy');
            cct = spd2cct(w,spd);
            il = illuminantSet(il,'name',sprintf('CCT %.0f',cct));
        end
        
        if isfield(illuminant,'comment'),
            il = illuminantSet(il,'comment',illuminant.comment);
        end
    end
end

end