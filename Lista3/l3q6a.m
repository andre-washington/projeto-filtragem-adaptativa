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

clear all
clf (figure(1))
clf (figure(2))

rep = 500; % numero de repeti��es
SNR = 30;
% n_var_4qam = 10^-3;

mu_n = 0.4;
% gama = 10^-9; 
   
% sinal 4 qam
constelation = [1+1i,1-1i,-1+1i,-1-1i];
data  = ceil(4*rand(rep,1));
s = transpose(constelation(data));  % sinal de treinamento

% n = n_var_4qam*randn(1,rep); %ru�do


H_num = [0.5 1.2 1.5 -1];
H_den = 1;
sh = filter(H_num, H_den, s);
n_var_4qam = var(sh) * 10^(-SNR/10);
    
%ru�do deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_4qam/2)*(randn(1,rep)+1i*(randn(1, rep)));
% transposte deve ser usada no lugar do hermitiano
x=sh+transpose(n); % sinal desejado na entrada do equalizador + ruido

% part I: trainning.
w(1,1:15) = 0; 
% xw(1,4) = 0;
xr = zeros(1, 15);
passo = 0.19;
gama = 10^(-6);
for k=15:rep
   
    xr = [x(k) xr(1:15-1)]; % varre os elementos de x
    xw = xr*w';
    e=s(k-14) - xw; %xr * w';
    erro(k)=e;
    w = w + passo * xr * conj(e) / (xr * xr' + gama);
end


figure(1)
plot(1:500,real(erro.*erro));
xlabel('plot(e.*e)');

% part II: decision block included. 

rep = 5000;
M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
dataIn = randi([0 1],k*rep,1);
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers

s2 = qammod(dataSymbolsIn,M,0); 

sh2 = filter(H_num, H_den, s2);
n_var_16qam = var(sh2) * 10^(-SNR/10);
    
%ru�do deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_16qam/2)*(randn(1,rep)+1i*(randn(1, rep)));
% transposte deve ser usada no lugar do hermitiano
x2=sh2+transpose(n); % sinal desejado na entrada do equalizador + ruido
 xr = zeros(1, 15);
 xd =  zeros(1,4995);
 for k=16:rep
    xr = [x2(k) xr(1:15-1)]; % varre os elementos de x
    xw = xr*w';
    
    % decision block
    xd(k) = round(real(xw)) + round(imag(xw))*1i ; 
   
    if(real(xd(k))>3)
        xd(k)=3+imag(xd(k))*1i;
    end
    if(real(xd(k))<-3)
        xd(k)= -3 + imag(xd(k))*1i;
    end
    if(imag(xd(k))>3)
        xd(k)= real(xd(k)) + 3i;
    end
    
    if(imag(xd(k))<-3)
        xd(k)= real(xd(k)) - 3i;
    end
    % end of decision block
    
    e= xd(k) - xw;
    erro(k)=e;
    w = w + passo * xr * conj(e) / (xr * xr' + gama);
    

 end
 
 figure(2)
plot(1:5000,real(erro.*erro),'red');
 axis([-1 rep -1 1]);
xlabel('plot(e.*e)');


[num, rate ]= biterr(abs(real(s2))',abs(real(xd)))
