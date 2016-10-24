[h13,h23,h123,w123] = SepCone(45,45,2,32);
y_out = zeros(360,90);


for azi = 1:360;
    for elv = 1:90;
        signal = siggen3d(1,azi,elv,32);
        [p_1, y_out(azi,elv)] = timefilt(signal,h13,h23);
    end
    azi
end