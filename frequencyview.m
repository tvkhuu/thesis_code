function [] = frequencyview(inpsig)
% Gets a frequency dependent 3d signal or filter
% shows the 2d contour plot over time of the 3rd dimension
% Assume its already fft shifted but still complex



% Taking the 3D FFT
inpsig=inpsig/max(max(max(abs(inpsig)))); % Normalizing the FFT


for k = 1:size(inpsig,3);
    hold all;
    contour(abs(inpsig(:,:,k)));
    grid on;
end


end

