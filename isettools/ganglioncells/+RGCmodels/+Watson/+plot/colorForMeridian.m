function color = colorForMeridian(meridianName)
    indexedMeridians = RGCmodels.Watson.constants.indexedMeridians;
    index = strcmp(indexedMeridians, meridianName);
    color = squeeze(RGCmodels.Watson.constants.meridianColors(index,:));
end

