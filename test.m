
Nx=32;
Nt=32;

wm=0.85*pi;
theta=60*pi/180;

epsilon=1*pi/180;

a=tan(theta-epsilon);
b=tan(theta+epsilon);

hn1n3=zeros(256,256);


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
        else
            fprintf('error\n');
        end
    end
end


H=fftshift(fftn(hn1n3));

H=H/max(max(abs(H)));

figure
mesh(abs(H'));view(2)


