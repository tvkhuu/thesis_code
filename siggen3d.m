function [ sig ] = siggen3d(amp,azi,elv,N )
% creates a 3d signal from a given azimuth and elevation.
% outputs as a 3d matrix (64,64,64)
% frequency of 120Hz

elevation = elv * pi/180;
azimuth = (azi) * pi/180;


n = 2*N+1;
r = linspace(-N,N,n);
%f = 20; %% Frequency of 120Hz
% n = 1:(length(r));
%fs = 1024;

xt = r .* sin(azimuth) .* cos(elevation);% .* 4*sin(2*pi*r*pi*f/fs);
yt = r .* cos(elevation) .* cos(azimuth);
t = r .* sin(elevation);
plot3(xt,yt,t);
axis([-N N -N N -N N]);


%%% Normalising values by rounding

xtn = round(xt)+(33*ones(1,length(xt)));
ytn = round(yt)+(33*ones(1,length(yt)));
tn = round(t)+(33*ones(1,length(t)));

sig = zeros(65,65,65);

for m = 1:length(r);
     sig(xtn(m),ytn(m),m) = 3; %turning it into a sinusoidal
end

