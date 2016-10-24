function [sigTime, sigFreq] = Siggen(azimuth, elevation, Nx, Ny, Ns, cen_freq, bandwidth,delay)
% Creates a 3d signal in both time and frequency domain. This is built for
% simulating a 32x32 plane with ts time samples

psi = elevation*pi/180; %elevation
phi = azimuth*pi/180; %azimuth
omega0= cen_freq*pi; %frequency cutoff
AA= bandwidth; % bandwidth of the signal
sigTime=zeros(Nx,Ny,Ns); %initialises the signal

% Computing the 3D received signal.
for n1=1:Nx % spatial domain x
    for n2=1:Ny % Spatial domain y
        for n3=1:Ns % Time domain
            %kk creates the line, w(n1,n2,n3) create a gaussian signal
            kk=sin(psi)*cos(phi)*n1 + sin(psi)*sin(phi)*n2 + n3-delay; % This is the 3D direction cosine argument
            sigTime(n1,n2,n3) = cos(omega0*kk)*exp(-AA*kk^2); % 3D signal %% 100 is the delay in time
        end
    end
end
sigFreq = fftshift(fftn(sigTime));
sigFreq = sigFreq/max(max(max(abs(sigFreq)))); % Normalisation

end

