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