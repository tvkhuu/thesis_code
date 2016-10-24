function [p1_ct ] = sigflow1(signal, filter)
% Applying the signal flow graph from Chamith's paper. Input signal is in
% 3d so doing this for each spatial and temporal while maintaining a
% contant single spatial output. This should output an array

sig_x = size(signal,1);
sig_y = size(signal,2);
sig_t = size(signal,3);
p1_ct = zeros(sig_y,sig_t);
y_ct = zeros(sig_t,1);
%filter is only in two dimensions

for t = 1:sig_t;
    for y = 1:sig_y;
        for Nct = 1:sig_t;
            x_sum = 0;
            for x = 1:sig_x;
                tmp = signal(x,y,t) .* filter(x,(sig_t+1)-Nct);
                x_sum = x_sum + tmp;
            end
            y_ct(Nct) = x_sum;
        end
        p1_ct(:,y,t) = y_ct;
    end
end

end

