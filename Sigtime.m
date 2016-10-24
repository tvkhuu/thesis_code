function [Sigt] = Sigtime(azimuth, elevation, ts, nt,d, bandwidth)
% Creates a 3d signal in both time and frequency domain. This is built for
% simulating a 32x32 plane with ts time samples

Nx=nt; % Number of antennas in the x-direction
Ny=nt; % Number of antennas in the y-direction
Nt=ts; % Number of sample points in the time

% Direction of arrivale in x-y-z spatial domain
psi= elevation*pi/180; % Angle from z-axis (elevation)   
phi= azimuth*pi/180; % Angle from x-axis (azimuth)

Sigt=zeros(Nt,Nt,Nt); % Storage place for 3D array signal. I have zero padded the signal in x and y dimensions to allow increased frequency resolution when taking FFT

% Parameters for the input signal. Here I used a cosine modulated Gaussian
% signal of the form x(t)=cos(omega0*t)exp(-AA*t^2), where omega0 sets the
% center frequency and AA sets the bandwidth. You may select any signal you
% want. The spectrum will have a line-shaped region of support and it is
% this signal that sets the actual spectral content on the line.
AA = bandwidth;
% Computing the 3D received signal.
for n1=1:Nx % spatial domain x
    for n2=1:Ny % Spatial domain y
        for n3=1:Nt % Time domain
            %kk creates the line, w(n1,n2,n3) create a sinusoindal signal
            kk= sin(psi)*cos(phi)*n1 + sin(psi)*sin(phi)*n2+n3-d; % This is the 3D direction cosine argument, n3 is the delay factor
            Sigt(n1,n2,n3) = exp(-AA*kk^2); % 3D signal impulse
        end
    end
end