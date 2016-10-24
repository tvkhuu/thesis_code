function [] = ContComp(signal, filter)
% Used to compare two filters on their contour intervals at level k
% Testing if the filter falls within the signal range
% done on the same figure

for k = 1:size(signal,3);
    clf
    contour(abs(signal(:,:,k)));
    hold on;
    contour(abs(filter(:,:,k)));
    grid on
    pause(0.10);
end


