
w0 = pi;
si=zeros(90,1);
phi=zeros(360,1);
X=zeros(90,360);
Y=zeros(90,360);
Z=zeros(90,360);
for aa=1:90
    %si(aa)=(aa-rr_si/2)*(pi/2)/(rr_si/2);
    si(aa)=(aa-90/2)*(pi/3)/(90/2);
    for bb=1:360
        phi(bb)=bb*(2*pi)/360;
        wx=w0*sin(si(aa))*cos(phi(bb));
        wy=w0*sin(si(aa))*sin(phi(bb));
        H(aa,bb)= arraypat(aa,bb);
        X(aa,bb)=abs(H(aa,bb))*sin(si(aa))*cos(phi(bb));
        Y(aa,bb)=abs(H(aa,bb))*sin(si(aa))*sin(phi(bb));
        Z(aa,bb)=abs(H(aa,bb))*cos(si(aa));
    end
end
%H=H/max(max(abs(H)));
mesh(X,Y,Z);