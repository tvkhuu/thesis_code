sig1 = siggen3d(3,45,45,32);
sig2 = siggen3d(2,20,65,32);
sig3 = siggen3d(1,90,20,32);

incoming = sig1 + sig2;

for t = 1:size(incoming,3);
    yn1 = conv2(incoming(:,:,t),hn1n3,'same');
    yn2(:,:,t) = conv2(yn1,hn2n3,'same');
end

ynoutput = abs(yn2); %% magnitude only

hold all
signalview(yn2);


%  yn1 = conv2(incoming(:,:,t),d3hn1n3(:,:,t),'same'); %% first fan filter
%  yn2(:,:,t) = conv2(y1,d3hn2n3(:,:,t),'same'); %% 2nd fan filter
