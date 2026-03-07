function defaultParams = coneMosaicTreeShrewDefaultParams()

    defaultParams = struct(...
        'centralRetinaConeSpacingMicrons', 6.0, ...   % Müller and Peichl, 1989
        'innerSegmentDiameterMicrons', 5.5, ...       % Müller and Peichl, 1989
        'coneDensities', [0.93 0 0.07] ...
    );

end
