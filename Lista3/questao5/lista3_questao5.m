
close
clear
%% Lista 3 - Questão 5
num_ite = 10000;
num_taps = 12;
mu_max = 1;
mu = [mu_max/2, mu_max/10, mu_max/50]; 
% input signal parameters
var_x = 1;
mean_x = 0;

x = zeros(num_ite,1);
%x = var_x*rand(num_coef,1) + mean_x;

% noise parameters
var_n = 10^-3;
mean_n = 0;

n = zeros(num_ite,1);
%n = var_n*randn(num_coef,1) + mean_n;

% filter taps
w = zeros(num_taps, 1);

%error vector
err = zeros(num_ite-12, 1);
bp=0;
y_sys = 0;
for i = 1:num_ite
    x(i) = var_x*rand + mean_x; 
    n(i) = var_n*randn + mean_n;
    if(i > 12)
        y_sys = x(i) - x(i-12) + y_sys;
        d = y_sys + n(i);
        
        y_flt = dot(x(i-11:i), w);
        err(i-12) = d - y_flt;
        w(i) = w(i-1) + 2*mu(2)*x(i-11:i)*err(i-12);
        if(err(i-12)^2 < 10^-8)
            bp = i;
            break;
        end
    end
    
end

plot(1:bp-12, err(1:bp-12).^2);