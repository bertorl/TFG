function [w, mse, e] = lms(N, mu, input, symbols, delay)

L = (N-1)/2;
q = 1;
w = zeros(N,1);
input = [zeros(1,L);input;zeros(1,L)];
    for k=1:length(symbols)
        ini=k-1+N-L;
        fin=k-1+N+L;
        samples = ini:1:fin;
        r_vector=input(samples);
        e(k)=symbols(k)-w'*r_vector;
        w=w+mu*e(k)*r_vector;
        mse = mean(e(k).^2);
    end
    mu = mu/k;
end