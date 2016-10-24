%%% array pattern
%%% cone is 60,30

Nx = 16; %number of X antennas
Ny = 16; %number of Y antennas
Nt = 16; %temporal Filter Order
Ns = 300; %number of time samples
ts = (1:Ns)';
delay = 40; %delay for original generated signal
psi = 30*pi/180; %elevation
phi = 60*pi/180; %azimuth
eps = 2*pi/180;  %fan angle
wm = pi; %Cut off frequency for filters
%omega0= 0.5*pi; %frequency cutoff
AA=0.05; % bandwidth of the signal
Wb = 16; %signal bit width
Fb = 12; %fractional bit position

elv = pi/2 - atan(sin(psi));
azi = phi;
theta_x = acos(cos(elv)*cos(azi));
theta_y = atan(tan(elv)/sin(azi));

a = tan(theta_x - eps);
b = tan(theta_x + eps);
c = tan(theta_y - eps);
d = tan(theta_y + eps);
for varia = 1:25
    omega0 = (varia*2+49)/100 * pi;
    
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
    outputs_elv_azi{90,360} = zeros(90,360);
    arraypat = zeros(90,360);
    for azi = 1:180;
        for elv = 1:90;
            sigTime = zeros(Nx,Ny,Ns);
            for n1=1:Nx % spatial domain x
                for n2=1:Ny % Spatial domain y
                    for n3=1:Ns % Time domain
                        %kk creates the line, w(n1,n2,n3) create a gaussian signal
                        kk=sin(elv*pi/180)*cos(azi*pi/180)*n1 + sin(elv*pi/180)*sin(azi*pi/180)*n2 + n3-delay; % This is the 3D direction cosine argument
                        sigTime(n1,n2,n3) = cos(omega0*kk)*exp(-AA*kk^2); % 3D signal %% 100 is the delay in time
                    end
                end
            end
            prefilt_sigti = cat(3,zeros(size(hn1n3,1), size(hn1n3,2), (size(hn1n3,2)-1)), sigTime);
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
            arraypat(elv,azi) = rms(time_output);
            elv
        end
        azi
    end
    arrfreq{varia} = arraypat;
end