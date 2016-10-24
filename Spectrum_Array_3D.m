
% Dr Chamith Wijenayake
% Date: 16/05/2016
% 3D space-time spectrum of a plane wave signal received by a planar array
% of Nx*Ny antennas

close all;
clc;

Nx=64; % Number of antennas in the x-direction
Ny=64; % Number of antennas in the y-direction
Nt=64; % Number of sample points in the time

% Direction of arrival in x-y-z spatial domain
psi=45*pi/180; % Angle from z-axis (elevation)   
phi=45*pi/180; % Angle from x-axis (azimuth)

w=zeros(Nt,Nt,Nt); % Storage place for 3D array signal. I have zero padded the signal in x and y dimensions to allow increased frequency resolution when taking FFT

% Parameters for the input signal. Here I used a cosine modulated Gaussian
% signal of the form x(t)=cos(omega0*t)exp(-AA*t^2), where omega0 sets the
% center frequency and AA sets the bandwidth. You may select any signal you
% want. The spectrum will have a line-shaped region of support and it is
% this signal that sets the actual spectral content on the line.
omega0=0.5*pi;
AA=0.1;

% Computing the 3D received signal.
for n1=1:Nx % spatial domain x
    for n2=1:Ny % Spatial domain y
        for n3=1:(Nt*4) % Time domain
            %kk creates the line, w(n1,n2,n3) create a sinusoindal signal
            kk=sin(psi)*cos(phi)*n1 + sin(psi)*sin(phi)*n2 + n3-200; % This is the 3D direction cosine argument
            w(n1,n2,n3) = cos(omega0*kk)*exp(-AA*kk^2); % 3D signal %% 100 is the delay in time
        end
    end
end

% Taking the 3D FFT
W=fftshift(fftn(w));
W=W/max(max(max(abs(W)))); % Normalizing the FFT

elv = pi/2 - atan(sin(psi));
azi = phi;
theta_x = acos(cos(elv)*cos(azi));
theta_y = atan(tan(elv)/sin(azi));
eps = 5 * pi/180; 

% move them to the 2 fan filter cases
a = tan(theta_x - eps);
b = tan(theta_x + eps);
c = tan(theta_y - eps);
d = tan(theta_y + eps);

x = linspace(-Nx/2,Nx/2,Nx+1);
y = linspace(-Ny/2,Ny/2,Ny+1);
z = linspace(-Nt/2,Nt/2,Nt+1);

%%% cut off frequency 
wm = pi;

%%% Create the two fan filters
hn1n3 = zeros(Nt,Nt);
hn2n3 = zeros(Nt,Nt);

for n1 = 1:Nx;
    for n3 = 1:Nt;
        if x(n1) == 0 && z(n3) ==0
            hn1n3(n1,n3) = ((b-a)*wm^2)/(4*a*b*pi^2);
        elseif x(n1) ~= 0 && ((x(n1)/a)+z(n3)) ~= 0 && ((x(n1)/b)+z(n3)) ~= 0
            hn1n3(n1,n3) = (1/(2*x(n1)*pi^2)) * ((1-cos(wm*((x(n1)/a)+z(n3))))/((x(n1)/a)+z(n3)) - ((1-cos(wm*((x(n1)/b)+z(n3))))/((x(n1)/b)+z(n3))));
        elseif ((x(n1)/b)+z(n3)) == 0 && ((x(n1)/a)+z(n3)) ~= 0;
            hn1n3(n1,n3) = (1/(2*x(n1)*pi^2))*((1-cos(wm*((x(n1)/a)+z(n3))))/((x(n1)/a)+z(n3)));
        elseif ((x(n1)/b)+z(n3)) ~= 0 && ((x(n1)/a)+z(n3)) == 0;
            hn1n3(n1,n3) = (1/(2*x(n1)*pi^2))*(-(1-cos(wm*((x(n1)/b)+z(n3))))/((x(n1)/b)+z(n3)));
        elseif x(n1) == 0 && z(n3) ~= 0
            hn1n3(n1,n3) = ((b-a)/(2*a*b*pi^2*z(n3)^2)) * (cos(wm*z(n3)) + wm*z(n3)*sin(wm*z(n3))-1);
        else
            fprintf('error\n');
        end
    end
end


for n2 = 1:Ny;
    for n3 = 1:Nt;
        if y(n2) == 0 && z(n3) ==0
            hn2n3(n2,n3) = ((d-c)*wm^2)/(4*c*d*pi^2);
        elseif y(n2) ~= 0 && ((y(n2)/c)+z(n3)) ~= 0 && ((y(n2)/d)+z(n3)) ~= 0
            hn2n3(n2,n3) = (1/(2*y(n2)*pi^2)) * ((1-cos(wm*((y(n2)/c)+z(n3))))/((y(n2)/c)+z(n3)) - ((1-cos(wm*((y(n2)/d)+z(n3))))/((y(n2)/d)+z(n3))));
        elseif ((y(n2)/d)+z(n3)) == 0 && ((y(n2)/c)+z(n3)) ~= 0;
            hn2n3(n2,n3) = (1/(2*y(n2)*pi^2))*((1-cos(wm*((y(n2)/c)+z(n3))))/((y(n2)/c)+z(n3)));
        elseif ((y(n2)/d)+z(n3)) ~= 0 && ((y(n2)/c)+z(n3)) == 0;
            hn2n3(n2,n3) = (1/(2*y(n2)*pi^2))*(-(1-cos(wm*((y(n2)/d)+z(n3))))/((y(n2)/d)+z(n3)));
        elseif y(n2) == 0 && z(n3) ~= 0
            hn2n3(n2,n3) = ((d-c)/(2*c*d*pi^2*z(n3)^2)) * (cos(wm*z(n3)) + wm*z(n3)*sin(wm*z(n3))-1);
        else
            fprintf('error\n');
        end
    end
end


% converting two time domain fans into the frequency domain
hw1w3 = fftshift(fft2(hn1n3));
hw2w3 = fftshift(fft2(hn2n3));


d3hw1w3 = zeros(Nt,Nt,Nt);
d3hw2w3 = zeros(Nt,Nt,Nt);

% extend both into each other's spatial domain to create wedge filters
for y = 1:size(hn2n3,1);
    d3hw1w3(:,y,:) = hw1w3(:,:);
end
for x = 1:size(hn1n3,1);
    d3hw2w3(x,:,:) = hw2w3(:,:);
end

% generate the 3d frequency overall domain through multiplying the two
% frequency wedges
Hw1w2w3 = d3hw2w3.*d3hw1w3;
%%% Normalising
Hw1w2w3 = Hw1w2w3/(max(max(max(abs(Hw1w2w3)))));
hw1w3 = hw1w3/(max(max(max(abs(hw1w3)))));
hw2w3 = hw2w3/(max(max(max(abs(hw2w3)))));
d3hw1w3 = d3hw1w3/(max(max(max(abs(d3hw1w3)))));
d3hw2w3 = d3hw2w3/(max(max(max(abs(d3hw2w3)))));


% figure
% plot(squeeze(w(Nx,Ny,:))); % 1D signal seen at the last antenna (i.e. Nx,Ny th antenna)
% grid on;
% title('1D time domain signal observed at a single (last) antenna')
% 
% figure
% mesh(abs(W(:,:,128+64))); % 3D spectrum visualized at some slice of omega3
% title('Magnitude frequency spectrum H(\omega_1,\omega_2,\omega_3=0.5\pi)')

%%% testing the filter with the overall cone
%sig_cone_out = W.*Hw1w2w3;

%%% Putting it through each fan cascaded
%sig_p_freq = W.*d3hw1w3;
%sig_cas_out = sig_p_freq.*d3hw2w3;

% Generating a signal with several DOAs
[sig1t,sig1w] = Siggen(75,15,Nt*4,0.5,0.1,100);
%[sig2t,sig2w] = Siggen(0, 90,Nt,0.4,0.1,30);
%[sig3t,sig3w] = Siggen(75,15,Nt,0.5,0.1,50);
sigsum_W = W+sig1w;
sigsum_T = w+sig1t;
%filt_cone_w = sigsum_W .* Hw1w2w3;

%filt_fan1_w = sigsum_W .* d3hw1w3;
%filt_out_w = filt_fan1_w .* d3hw2w3;


%Do the time domain implementation

%Visualizing the 3D line shaped region of support
w1=((-pi+2*pi/Nt):2*pi/Nt:pi);
w2=((-pi+2*pi/Nt):2*pi/Nt:pi);
w3=((-pi+2*pi/(Nt*4)):2*pi/(Nt*4):pi);
[W1,W2,W3] = meshgrid(w1,w2,w3);
figure
p = patch(isosurface(W1,W2,W3, abs(sigsum_W),0.7));
isonormals(W1,W2,W3, abs(sigsum_W), p)
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');
%daspect([1 1 1])
axis([-pi pi -pi pi -pi pi])
view(-47,12);
grid on
camlight; lighting phong
title('3D spectrum with Cone Filter')
hold on
f = patch(isosurface(W1,W2,W3, abs(Hw1w2w3), 0.5));
isonormals(W1,W2,W3, abs(Hw1w2w3), f);

[p_1, y_out] = timefilt(sigsum_T,hn1n3,hn2n3);
figure
mesh(p_1);
title('After one Fan filter');
xlabel('Spatial Position');
ylabel('Samples over time');

figure
plot(y_out);
title('After both filters');
xlabel('Time per sample');


figure
mesh(abs(sigsum_W(:,:,192)));