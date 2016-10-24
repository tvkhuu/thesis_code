function [p_t, time_output] = timefilt(signal, filter_1, filter_2)
% Time Domain filtering using a signal flow graph. Cascading both fan
% filters
% input signal has three dimensions
% Filter across one dimension so y constant
% padding the signal by zeros
% do this for one point in time

prefilt_sigti = cat(3,zeros(size(filter_1,1), size(filter_1,2), (size(filter_1,2)-1)), signal);
p_t = zeros(size(filter_1,1), size(filter_1,2));
for t = size(filter_1,2):size(prefilt_sigti,3);
    for y = 1:size(prefilt_sigti,2);
        h_1 = 0;
        for Nct = 1:size(filter_1,2);
            for x = 1:size(prefilt_sigti,1);
                tmp = prefilt_sigti(x,y,(t-(Nct-1))) * filter_1(x,size(filter_1,2)-(Nct-1));
                h_1 = h_1 + tmp;
            end
        end
        p_t(y,t-(size(filter_1,2)-1)) = h_1;
    end
end
time_output = zeros(size(p_t,2),1);
prefilt_p_t = cat(2,zeros(size(filter_2,1), size(filter_2,2)-1), p_t);
for t = size(filter_2,2):size(prefilt_p_t,2);
    h_2 = 0;
    for Nct = 1:size(filter_2,2);
        for y = 1:size(filter_2,1);
            tmp = prefilt_p_t(y,(t-(Nct-1))) * filter_2(y,size(filter_2,2)-(Nct-1));
            h_2 = h_2 + tmp;
        end
    end
    time_output(t-(size(filter_2,2)-1)) = h_2;
end