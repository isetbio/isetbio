function defaultParams = coneMosaicTreeShrewDefaultParams()

    defaultParams = struct(...
        'centralRetinaConeSpacingMicrons', 6.0, ...   % Müller and Peichl, 1989
        'innerSegmentDiameterMicrons', 5.5, ...       % Müller and Peichl, 1989
        'outerSegmentLengthMicrons', 10, ...          % SAMORAJSKI , ORDY , and KEEFE, Structural organization of the retina in the tree shre (tupaia glis), THE JOURNAL OF CELL BIOLOGY • VOLUME ~8, 1966
        'coneDensities', [0.93 0 0.07] ...
    );

end
