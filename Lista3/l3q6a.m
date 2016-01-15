% % LMS complexo
% e(k) = d(k) -w^H (k) x(k);
% w(k+1) = w(k) + mu_c * conjugated(e) * (k) * x(k)

% % NRMS
% % Inicializa
% x(0) = w(0) = [ 0 0 .. 0]^T
% mu_n no intervalo 0< mu_n <= 2
% gama  const. pequeno
% 
% % para k> 0:
% e(k) = d(k) - x^T (k) *w(k)
% w(k+1) = w(k) + mu_n/ (gama + x^T(k) * x(k) ) * e(k) * x(k)

clear
rep = 500; % numero de repetições
SNR = 30;
% n_var_4qam = 10^-3;

mu_n = 1;
gama = 10^-9; 
   
% sinal 4 qam
const = [1+i,1-i,-1+i,-1-i];
data  = ceil(4*rand(rep,1));
s = transpose(const(data));  % sinal de treinamento

% n = n_var_4qam*randn(1,rep); %ruído


H_num = [0.5 1.2 1.5 -1];
H_den = 1;
sh = filter(H_num, H_den, s);
n_var_4qam = var(sh) * 10^(-SNR/10);
    
%ruído deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_4qam/2)*(randn(1,rep)+1i*(randn(1, rep)));
% transposte deve ser usada no lugar do hermitiano
x=sh+transpose(n); % sinal desejado na entrada do equalizador + ruido

w(1,1:4) = 0; 
xw(1,4) = 0;
xr = zeros(1, 4);
passo = 0.1;
gama = 10^(-6);
for i=4:rep
    xr = [x(i) xr(1:4-1)]; % varre os elementos de x
    s(i-3);
    e=s(i-3) - xr * w';
    w = w + passo * xr * conj(e) / (xr * xr' + gama);
   

hold on
semilogx(i,real(w));
xlabel('plot(w)');
end
% subplot(2,1, 1);
% semilogx(real(e.*e));
% xlabel('plot(e.*e)');
% 
% subplot(2,1, 2);
% plot(real(e));
% xlabel('plot(e)');