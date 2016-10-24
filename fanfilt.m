%%% creating a fan filter to generate a generic response to compare it to
%%% the inputs of the verilog simulation to see if it is similar

Nx=8;
Nt=8;

wm=pi;

theta=45*pi/180;
epsilon=10*pi/180;

a=tan(theta-epsilon);
b=tan(theta+epsilon);
   
fhn1n3=zeros(Nx,Nt);

for n1 = 1:Nx;
    x=(n1-Nx/2);
    for n3 = 1:Nt;
         z=(n3-Nt/2);
        if x == 0 && z ==0
            fhn1n3(n1,n3) = ((b-a)*wm^2)/(4*a*b*pi^2);
        elseif x ~= 0 && ((x/a)+z) ~= 0 && ((x/b)+z) ~= 0
            fhn1n3(n1,n3) = (1/(2*x*pi^2)) * ((1-cos(wm*((x/a)+z)))/((x/a)+z) - ((1-cos(wm*((x/b)+z)))/((x/b)+z)));
        elseif ((x/b)+z) == 0 && ((x/a)+z) ~= 0;
            fhn1n3(n1,n3) = (1/(2*x*pi^2))*((1-cos(wm*((x/a)+z)))/((x/a)+z));
        elseif ((x/b)+z) ~= 0 && ((x/a)+z) == 0;
            fhn1n3(n1,n3) = (1/(2*x*pi^2))*(-(1-cos(wm*((x/b)+z)))/((x/b)+z));
        elseif x == 0 && z ~= 0
            fhn1n3(n1,n3) = ((b-a)/(2*a*b*pi^2*z^2)) * (cos(wm*z) + wm*z*sin(wm*z)-1);
        else
            fprintf('error\n');
        end
    end
end

fhw1w3 = fftshift(fftn(fhn1n3));
% normalising
fhw1w3 = fhw1w3/max(max(abs(fhw1w3)));
fhn1n3 = fhn1n3/max(max(abs(fhn1n3)));

% signal flow graph
p_1t = zeros(Nt,1);

%generate 2d signal using a slice from the 3d signal
% implement the signal flow graph and see the output as the vector passing
% at the filter DOA
for t = 1:size(sig2dt,2);
    h_1 = zeros(size(fhn1n3,2),1);
    for Nct = 1:size(fhn1n3,2);
        for x = 1:size(sig2dt,1);
            tmp = sig2dt(x,t) * fhn1n3(x,size(fhn1n3,2)-Nct+1);
            h_1(Nct) = h_1(Nct) + tmp; 
        end
    end
    p_1t(t) = sum(h_1);
end

figure
plot(p_1t);



