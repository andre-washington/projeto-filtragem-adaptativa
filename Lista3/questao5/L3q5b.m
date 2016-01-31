% versao 2
close
clear

rep = 1000; % numero de repetições
% mu_max = 0.001;

x_var = 1; %variancia de x
n_var = 10^-3;% variancia de n

x = x_var*rand(1,rep)+5; %sinal de entrada
n = n_var*randn(1,rep); %ruído

w(1:12,1) =0;
h_num = [1 0 0 0 0 0 0 0 0 0 0 0 -1];
h_den = [1 -1];
yh = filter(h_num, h_den, x);

mu_max = 1/trace(x_var.*eye(12));
mu = [mu_max/2, mu_max/10, mu_max/50]; 

for i = length(w) : rep
    x_v = fliplr(yh(i-11:i))'
    d = yh(i-11)+n(i-11);
    e = d - w'*x_v;
    
    w = w + 2 * mu_max * e * x_v;

    mse(i-11)  = e.*e;
 end

figure(1)

semilogy(mse);
title('LMS'); xlabel('amostras'); ylabel('MSE');
