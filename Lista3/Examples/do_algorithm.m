function [mse] = do_algorithm (samples, data, x , w, p, R, alg)
mse = zeros(1, length(samples));
w_o=w;
for i=length(w):samples
    x_v = [x(i); x(i-1)];
    d = data(i-1);
    e   = d - w'*x_v;
    
    switch lower(alg)
        case {'lms'}
            w = w + 2 * 0.01 * e * x_v;
        case {'newton'}
            w = w - 0.1 * (R \ (-2 * p + 2 * R * w));
        case {'nlms'}
            w = w + 0.02 * (1/(x_v' * x_v + 0.01)) * e.* x_v;
        case {'grad'}
            w = w - 2 * 0.01 * (R * w - p);
        %case{'rls'}
    end
    plot([w_o(1) w(1)], [w_o(2) w(2)]);
    w_o = w;
    mse(i-1)  = e.*e;
    
end