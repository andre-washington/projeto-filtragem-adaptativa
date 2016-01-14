% clc
clear

rep = 1000; % numero de repetições
mu_max = 0.001;
x_var = 1; %variancia de x
n_var = 10^-3;% variancia de n
x = x_var*randn(1,rep)+5; %sinal de entrada
n = n_var*randn(1,rep); %ruído
w(1:12,1) =0;

yh(1) = x(1);


for i=1:rep
    if(i>12)
        yh(i) = x(i) - x(i-12) + yh(i-1); % saida do sistema H
    
    else
        yh(i+1) = x(i+1) - 0 + yh(i);

    end

    d(i) = yh(i) + n(i);
    
    
    if(i>12)
%         dot(x(i-11:i), w);
        yw(i) = dot(w(1:12),x(i-(1:12)));
        e(i) = d(i) - yw(i);
        w(1:12) =  w(1:12)' + 2*mu_max*e(i).* x(i-(1:12)) ;
    else
        yw(i) = sum(w(1:i)) ;
        e(i) = d(i) - yw(i);
        w(1:12) = w(1:12) + 2*mu_max*e(i);
        

    end
%     figure(1);
%     hold on
%     semilogx(i,w);
    
    
end
% figure(2)

subplot(2,1, 1);

semilogy(e.*e);
xlabel('plot(e.*e)');
subplot(2,1, 2);
plot(e);
xlabel('plot(e)');




