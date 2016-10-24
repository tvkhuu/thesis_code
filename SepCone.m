function [hn1n3,hn2n3] = SepCone(azimuth,elevation,eps,Nx, Nt, wm)
% % From the given three angles, create a dual-spatial temporal cone filter
% % with angles, azimuth, elevation and cone angle by using a seperable
% % approach of two wedge filters. The two wedge filters are cascaded to
% % create the cone filter. 
% % Inputs Azimith, Elevation and cone width with number of antennas
% % Angles in degrees 
% % t is the order of the filter, number of time samples
% % Outputs Two wedge filters and the cone filter in discrete time.
% %%% Tran Vinh Khuu z3418060 Thesis B 2016 %%%
% the cone will have both x and y components which are multiplied these x
% and y components can be represented as a 2d fan filter on the plane or a
% wedge in 3d
% From a cone's azimuth and elevation we obtain the two x and y component
% angles

elv = pi/2 - atan(sin(elevation));
azi = azimuth;
theta_x = acos(cos(elv)*cos(azi));
theta_y = atan(tan(elv)/sin(azi));

a = tan(theta_x - eps);
b = tan(theta_x + eps);
c = tan(theta_y - eps);
d = tan(theta_y + eps);

hn1n3 = zeros(Nx,Nt);
hn2n3 = zeros(Nx,Nt);

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

for n2 = 1:Nx;
    y=(n2-Nx/2);
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
hn1n3 = hn1n3/max(max(abs(hn1n3)));
hn2n3 = hn2n3/max(max(abs(hn2n3)));

% hw1w3 = fftshift(fft2(hn1n3));
% hw2w3 = fftshift(fft2(hn2n3));
% hw1w3 = hw1w3/max(max(abs(hw1w3)));
% hw2w3 = hw2w3/max(max(abs(hw2w3)));


end
