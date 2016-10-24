%close all
%Parameters
Nx = 64; %number of X antennas
Ny = 64; %number of Y antennas
Nt = 64; %temporal Filter Order
Ns = 64; %number of time samples
ts = (1:Ns)';
delay = 50; %delay for original generated signal
psi = 30*pi/180; %elevation
phi = 60*pi/180; %azimuth
eps = 5*pi/180;  %fan angle
wm = pi; %Cut off frequency for filters
omega0=0.5*pi; %frequency cutoff
AA=0.05; % bandwidth of the signal

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

% signal generation or creating the signal
sigTime=zeros(Nt,Nt,Ns); %initialises the signal



for pv = 1:100;
    omega0 = (pv/100)*pi;
    for qv = 1:360;
        phi = qv*pi/180;
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
        
        hw1w3 = fftshift(fft2(hn1n3));
        hw2w3 = fftshift(fft2(hn2n3));
        hw1w3 = hw1w3/max(max(abs(hw1w3)));
        hw2w3 = hw2w3/max(max(abs(hw2w3)));
        d3hw1w3 = zeros(Nx,Ny,Nt);
        d3hw2w3 = zeros(Nx,Ny,Nt);
        for y = 1:size(hn2n3,1);
            d3hw1w3(:,y,:) = hw1w3(:,:);
        end
        for x = 1:size(hn1n3,1);
            d3hw2w3(x,:,:) = hw2w3(:,:);
        end
        hw1w2w3 = d3hw1w3.*d3hw2w3;
        filteroutfreq = hw1w2w3.*sigFreq;
        testfreq(qv,pv) = rms(filteroutfreq))));
        qv
    end
    pv
end


