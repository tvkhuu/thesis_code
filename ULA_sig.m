function [w] = ULA_sig(psi,N,x)
%The function accepts 3 arguments:
% 1. psi -- the angle of arrival in degree
% 2. N  --- number of antennas in the array
% 3. x ---- the interference signal
% return the 2D discrete representation of the interference w(n1,n2)
    Tau  = sind(psi); % fixed time delay between consecutive elements
    w(:,1) = x;
    for p = 2:N
        w(:,p) = FFT_delay(x,(p-1)*Tau);
    end
end

