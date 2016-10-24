function [ ] = signalview(sig)
% views a signal with time dependency as the third dimension
% meshes each slice of a 3d array over time
clf

sig =sig/max(max(max(abs(sig)))); % Normalizing the signal


for t = 1:size(sig,3);
    mesh(abs(sig(:,:,t)));
    colormap jet;
    grid on;
    axis([0 size(sig,1) 0 size(sig,2) 0 1]);
    view(10,80);
    pause(0.05);
end