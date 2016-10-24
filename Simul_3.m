close all
clear
%Parameters
Nx = 16; %number of X antennas
Ny = 16; %number of Y antennas
Nt = 16; %temporal Filter Order
Ns = 128; %number of time samples
ts = (1:Ns)';   
delay = 50; %delay for original generated signal
psi = 30*pi/180; %elevation
phi = 60*pi/180; %azimuth
eps = 5*pi/180;  %fan angle
wm = pi; %Cut off frequency for filters
omega0=0.5*pi; %frequency cutoff
AA=0.05; % bandwidth of the signal
intpsi = 85*pi/180; %interference elevation
intphi = 5*pi/180; %interference azimuth
intde = 100; %delay for interference signal
intomeg = 0.5*pi; %frequency cutoff for interference
intbw = 0.03; %bandwidth interference
Wb = 16; %signal bit width
Fb = 14; %fractional bit position

% signal generation or creating the signal
sigTime=zeros(Nt,Nt,Ns); %initialises the signal

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

%creating the filter
elv = pi/2 - atan(sin(psi));
azi = phi;
theta_x = acos(cos(elv)*cos(azi));
theta_y = atan(tan(elv)/sin(azi));

a = tan(theta_x - eps);
b = tan(theta_x + eps);
c = tan(theta_y - eps);
d = tan(theta_y + eps);

hn1n3 = zeros(Nx,Nt);
hn2n3 = zeros(Ny,Nt);

for n1 = 1:Nx;
    x=(n1-Nx/2);
    for n3 = 1:Nt;
        z=(n3-Nt/2);
        if x == 0 && z ==0
            hn1n3(n1,n3) = ((b-a)*wm^2)/(4*a*b*pi^2);
        elseif x ~= 0 && ((x/a)+z) ~= 0 && ((x/b)+z) ~= 0
            hn1n3(n1,n3) = (1/(2*x*pi^2)) * ((1-cos(wm*((x/a)+z)))/((x/a)+z) - ((1-cos(wm*((x/b)+z)))/((x/b)+z)));
        elseif ((x/b)+z) == 0 && ((x/a)+z) ~= 0;
            hn1n3(n1,n3) = (1/(2*x*pi^2))*((1-cos(wm*((x/a)+z)))/((x/a)+z));
        elseif ((x/b)+z) ~= 0 && ((x/a)+z) == 0;
            hn1n3(n1,n3) = (1/(2*x*pi^2))*(-(1-cos(wm*((x/b)+z)))/((x/b)+z));
        elseif x == 0 && z ~= 0
            hn1n3(n1,n3) = ((b-a)/(2*a*b*pi^2*z^2)) * (cos(wm*z) + wm*z*sin(wm*z)-1);
        end
    end
end

for n2 = 1:Ny;
    y=(n2-Ny/2);
    for n3 = 1:Nt;
        z=(n3-Nt/2);
        if y == 0 && z ==0
            hn2n3(n2,n3) = ((d-c)*wm^2)/(4*c*d*pi^2);
        elseif y ~= 0 && ((y/c)+z) ~= 0 && ((y/d)+z) ~= 0
            hn2n3(n2,n3) = (1/(2*y*pi^2)) * ((1-cos(wm*((y/c)+z)))/((y/c)+z) - ((1-cos(wm*((y/d)+z)))/((y/d)+z)));
        elseif ((y/d)+z) == 0 && ((y/c)+z) ~= 0;
            hn2n3(n2,n3) = (1/(2*y*pi^2))*((1-cos(wm*((y/c)+z)))/((y/c)+z));
        elseif ((y/d)+z) ~= 0 && ((y/c)+z) == 0;
            hn2n3(n2,n3) = (1/(2*y*pi^2))*(-(1-cos(wm*((y/d)+z)))/((y/d)+z));
        elseif y == 0 && z ~= 0
            hn2n3(n2,n3) = ((d-c)/(2*c*d*pi^2*z^2)) * (cos(wm*z) + wm*z*sin(wm*z)-1);
        end
    end
end
%Normailising
hw1w3 = fftshift(fft2(hn1n3));
hw2w3 = fftshift(fft2(hn2n3));
hw1w3 = hw1w3/max(max(abs(hw1w3)));
hw2w3 = hw2w3/max(max(abs(hw2w3)));
figure
mesh(abs(hw1w3));
title('Fan Filter Hn1n3 in frequency domain');

d3hw1w3 = zeros(Nx,Ny,Nt);
d3hw2w3 = zeros(Nx,Ny,Nt);
for y = 1:size(hn2n3,1);
    d3hw1w3(:,y,:) = hw1w3(:,:);
end
for x = 1:size(hn1n3,1);
    d3hw2w3(x,:,:) = hw2w3(:,:);
end
hw1w2w3 = d3hw1w3.*d3hw2w3;

w1=((-pi+2*pi/Nx):2*pi/Nx:pi);
w2=((-pi+2*pi/Ny):2*pi/Ny:pi);
w3=((-pi+2*pi/(Ns)):2*pi/(Ns):pi);
[W1,W2,W3] = meshgrid(w1,w2,w3);
figure
p = patch(isosurface(W1,W2,W3, abs(sigFreq),0.7));
isonormals(W1,W2,W3, abs(sigFreq), p)
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
%daspect([1 1 1])
axis([-pi pi -pi pi -pi pi])
view(-47,12);
grid on
camlight; lighting phong
title('3D spectrum with Cone Filter')
hold on
w1=((-pi+2*pi/Nx):2*pi/Nx:pi);
w2=((-pi+2*pi/Ny):2*pi/Ny:pi);
w3=((-pi+2*pi/(Nt)):2*pi/(Nt):pi);
[W1,W2,W3] = meshgrid(w1,w2,w3);
f = patch(isosurface(W1,W2,W3, (abs(hw1w2w3).^2), 0.3));
isonormals(W1,W2,W3, abs(hw1w2w3), f);
set(f, 'FaceColor', 'red', 'EdgeColor', 'none');

figure
contour(abs(hw1w2w3(:,:,(3*Nt/4))));
hold on
contour(abs(sigFreq(:,:,(3*Ns/4))));
title('Contour Slice at pi/2 matching filter to signal DoA in Frequency Domain');

%Adding interference to see if filter works
timeinter = zeros(Nt,Nt,Ns);
for n1=1:Nx % spatial domain x
    for n2=1:Ny % Spatial domain y
        for n3=1:Ns % Time domain
            %kk creates the line, w(n1,n2,n3) create a gaussian signal
            kk=sin(intpsi)*cos(intphi)*n1 + sin(intpsi)*sin(intphi)*n2 + n3-intde; % This is the 3D direction cosine argument
            timeinter(n1,n2,n3) = cos(intomeg*kk)*exp(-intbw*kk^2); % 3D signal %% 100 is the delay in time
        end
    end
end
freqinter = fftshift(fftn(timeinter));
freqinter = freqinter/max(max(max(abs(freqinter)))); % Normalisation

%Combining both
sigwintti = sigTime + timeinter;
sigwintfr = freqinter + sigFreq;

%Interference signal as seen from last antenna
figure
plot(squeeze(sigwintti(Nx,Ny,:)));
title('Signal seen at last antenna');
figure
mesh(abs(sigwintfr(:,:,(3*Ns/4))));
title('Interference and Desired signal at slice pi/2');

%filtering the interference out
% freq_output = hw1w2w3.*sigwintfr;
% figure
% mesh(abs(freq_output(:,:,(3*Ns/4))));
% title('Frequency filtering output');

%signal flow graph
%padding the signal with zero condition zeros by length of temportal order filter - 1
prefilt_sigti = cat(3,zeros(size(hn1n3,1), size(hn1n3,2), (size(hn1n3,2)-1)), sigwintti);

p_t = zeros(size(hn1n3,1), size(hn1n3,2));
for t = size(hn1n3,2):size(prefilt_sigti,3);
    for y = 1:size(prefilt_sigti,2);
        h_1 = 0;
        for Nct = 1:size(hn1n3,2);
            for x = 1:size(prefilt_sigti,1);
                tmp = prefilt_sigti(x,y,(t-(Nct-1))) * hn1n3(x,size(hn1n3,2)-(Nct-1));
                h_1 = h_1 + tmp;
            end
        end
        p_t(y,t-(size(hn1n3,2)-1)) = h_1;
    end
end
time_output = zeros(size(p_t,2),1);
prefilt_p_t = cat(2,zeros(size(hn2n3,1), size(hn2n3,2)-1), p_t);
for t = size(hn2n3,2):size(prefilt_p_t,2);
    h_2 = 0;
    for Nct = 1:size(hn2n3,2);
        for y = 1:size(hn2n3,1);
            tmp = prefilt_p_t(y,(t-(Nct-1))) * hn2n3(y,size(hn2n3,2)-(Nct-1));
            h_2 = h_2 + tmp;
        end
    end
    time_output(t-(size(hn2n3,2)-1)) = h_2;
end

figure
mesh(p_t);
title('Filtered first stage signal in time domain');
figure
mesh(abs(fftshift(fftn(p_t))));
title('Frequency Mesh of first stage signal from time domain');
figure
plot(time_output)
title('Filtered signal as seen from an antenna');

%%% Adjustments and allocations for verilog implementation
%sig_input = sigwintti;
sig_input = squeeze(sigwintti(:,4,:))';