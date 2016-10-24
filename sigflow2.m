function [y_out] = sigflow2(signal1, filter_2)
% The second filter to be cascaded to create the cone filter implementation
% as per Chamith's paper. The signal input is the output from the function
% sigflow1 and the filter being the 2nd fan filter

sig_x = size(signal1,1);
sig_y = size(signal1,2);
sig_t = size(signal1,3);
y_ct = zeros(sig_t,1);
y_out = zeros(sig_y,sig_t);

for t = 1:sig_t;
    for Nct = 1:sig_t;
        sum_y = 0;
        for y = 1:sig_x;
            tmp = signal1(y,t) .* filter_2(y,(sig_t+1)-Nct);
            sum_y = sum_y + tmp;
        end
        y_ct(Nct) = sum_y;
    end
    y_out(:,t) = y_ct;
end

