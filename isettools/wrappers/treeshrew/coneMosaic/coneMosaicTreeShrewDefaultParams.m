function defaultParams = coneMosaicTreeShrewDefaultParams()

    defaultParams = struct(...
        'centralRetinaConeSpacingMicrons', 6.5, ...  % Müller and Peichl, 1989
        'innerSegmentDiameterMicrons', 6.0, ...       % Müller and Peichl, 1989
        'coneDensities', [0.9 0 0.1] ...
    );

end
