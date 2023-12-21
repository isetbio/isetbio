function fDir = figureDir(theParentFilename, saveFigures)
    fDir = fullfile(isetbioRootPath,'local',theParentFilename);
    if (saveFigures)
        if (~exist(fDir,'dir'))
            mkdir(fDir);
        end
        fprintf('Will save figures into %s\n',fDir);
    else
        fprintf('Not saving figures. Set saveFigures to true in ''%s.m'' to save\n', theParentFilename);
    end
end