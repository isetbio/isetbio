function visualizeDisplaySPD(wavelengthSupport, spds)
% Method to visualize the display's SPDs
figure(); clf;
% plot the Red primary SPD
bar(wavelengthSupport, spds(:,1), 1, 'FaceColor', [1 0.5 0.5], 'EdgeColor', [1 0 0], 'FaceAlpha', 0.5); hold on
% plot the Green and Blue primary SPDs
bar(wavelengthSupport, spds(:,2), 1, 'FaceColor', [0.5 1 0.5], 'EdgeColor', [0 1 0], 'FaceAlpha', 0.5);
bar(wavelengthSupport, spds(:,3), 1, 'FaceColor', [0.5 0.5 1], 'EdgeColor', [0 0 1], 'FaceAlpha', 0.5);
% label plot
xlabel('\lambda (nm)'); 
ylabel('\it radiance (Watts  \cdot sr^{-1} \cdot m^{-2} \cdot nm^{-1})'); 
set(gca, 'FontSize', 16);
title('spectral radiance')
legend({'R primary', 'G primary', 'B primary'});
end
