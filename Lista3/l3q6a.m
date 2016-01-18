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
%
% para o caso complexo:
% e(k) = d(k) - x^H (k) *w(k)
% w(k+1) = w(k) + mu_n/ (gama + x^H(k) * x(k) ) * conjugated (e(k)) * x(k)


clear all
clf (figure(1))
clf (figure(2))

num_taps = 15;
rep = 500; % numero de repetições
SNR = 30;
% n_var_4qam = 10^-3;

mu = 0.4;
gama = 10^-9; 
   
% sinal 4 qam
%constelation = [1+1i,1-1i,-1+1i,-1-1i];
%data  = ceil(4*rand(rep,1));
%s = transpose(constelation(data));  % sinal de treinamento

M = 4;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol

rng default 
dataIn = randi([0 1],k*rep,1);
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers
s = qammod(dataSymbolsIn,M,0); 
% n = n_var_4qam*randn(1,rep); %ruído


H_num = [0.5 1.2 1.5 -1];
H_den = 1;
sh = filter(H_num, H_den, s);
n_var_4qam = var(sh) * 10^(-SNR/10);
    
%ruído deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_4qam/2)*(randn(rep, 1)+1i*(randn(rep, 1)));
% transposte deve ser usada no lugar do hermitiano
x=sh+n; % sinal desejado na entrada do equalizador + ruido

% part I: trainning.
w = zeros(num_taps, 1); 
% xw(1,4) = 0;
%xr = zeros(1, 4);

eq_out = zeros(1, rep);
err_vec = zeros(1,rep);

for k = num_taps:rep
    input = x(k-num_taps+1:k);
    eq_out(k) = w'*input;
    err_vec(k) = s(k) - eq_out(k);
    w = w + mu * conj(err_vec(k)) * input / ( input' * input + gama);
    
end


figure(1)
semilogy(1:500, real(conj(err_vec).*err_vec));
xlabel('plot(e.*e)');

% part II: decision block included. 

rep = 5000;
% sinal 16 qam
%constelation2 = [
%    1+1i,1-1i,-1+1i,-1-1i, 2+1i,2-1i,-2+1i,-2-1i, 1+2i,1-2i,-1+2i,-1-2i, 2+2i,2-2i,-2+2i,-2-2i ];
%data  = ceil(16*rand(rep,1));
%s2 = transpose(constelation2(data));  % sinal de entrada

M = 16;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
dataIn = randi([0 1],k*rep,1);
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers

s2 = qammod(dataSymbolsIn,M,0); 

sh2 = filter(H_num, H_den, s2);
n_var_4qam = var(sh2) * 10^(-SNR/10);
    
%ruído deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_4qam/2)*(randn(rep, 1)+1i*(randn(rep, 1)));
% transposte deve ser usada no lugar do hermitiano
x2=sh2+n; % sinal desejado na entrada do equalizador + ruido
 
eq_out = filter(w, 1, x2);
dataSymbolsOut = qamdemod(eq_out,M); 
dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:); 

[numErrors,ber] = biterr(dataIn,dataOut);

ber

% figure(2)
%plot(1:5000,real(erro.*erro),'red');
% axis([-1 600 -1 1]);
%xlabel('plot(e.*e)');

