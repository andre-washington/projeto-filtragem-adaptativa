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

clear

num_taps = 15;
rep = 500; % numero de repetiï¿½ï¿½es
SNR = 30;


mu = 0.4;
gama = 10^-5; 
   
% sinal 4 qam
%constelation = [1+1i,1-1i,-1+1i,-1-1i];
%data  = ceil(4*rand(rep,1));
%s = transpose(constelation(data));  % sinal de treinamento

M = 4;                     % Size of signal constellation
k = log2(M);                % Number of bits per symbol
 
dataIn = randi([0 1],k*rep,1);
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);   % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                 % Convert to integers
s = qammod(dataSymbolsIn,M,0); 
% n = n_var_4qam*randn(1,rep); %ruï¿½do


H_num = [0.5 1.2 1.5 -1];
H_den = 1;
sh = filter(H_num, H_den, s);
n_var_4qam = var(sh) * 10^(-SNR/10);
    
%ruï¿½do deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_4qam/2)*(randn(rep, 1)+1i*(randn(rep, 1)));
% transposte deve ser usada no lugar do hermitiano
x=sh+n; % sinal desejado na entrada do equalizador + ruido

% part I: trainning.
w = zeros(num_taps, 1); 

eq_out = zeros(1, rep);
err_vec = zeros(1,rep);

init = zeros(num_taps - 1, 1);

% for k = 1:rep
%     if(k < num_taps)
% 	inp = [x(k); init];
% 	init = inp(1:end-1);
%     else
% 	inp = x(k:-1:k-num_taps+1);	   
%     end
%     eq_out(k) = w'*inp;
%     err_vec(k) = s(k) - eq_out(k);
%     w = w + (mu/(gama + inp' * inp))*(conj(err_vec(k)) * inp);
%     
% end

for k = 8:rep % Ã‰ nessessario esperar um momento atÃ© fazer a comparaÃ§Ã£o com o sinal (SeÃ§Ã£o 2.10.4 pag.57 Diniz), uma boa espera Ã© a metade do comprimento do filtro. Nesse caso 7 ou 8.
    if(k < num_taps)
	inp = [x(k); init];
    init = inp(1:end-1);
    else
	inp = x(k:-1:k-num_taps+1);	   
    end
    eq_out(k) = w'*inp;
    err_vec(k) = s(k-7) - eq_out(k); % Compara o sinal de entrada com saida de saida (adiantado). Por exemplo(sinal 1 com saÃ­da 8)
    w = w + (mu/(gama + inp' * inp))*(conj(err_vec(k)) * inp);
    
end
figure(1)
plot3(real(eq_out),imag(eq_out),1:rep,'r.'); % plota da modulaÃ§Ã£o 4 QAM
%keyboard;
figure(2)
semilogy(1:500, real(conj(err_vec).*err_vec));
xlabel('plot(e.*e)');

%figure(2)
%semilogy(real(s));
%xlabel('plot(e.*e)');



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
    
%ruï¿½do deve ser complexo e a variancia divida para cada parte (real/imag)
n = sqrt(n_var_4qam/2)*(randn(rep, 1)+1i*(randn(rep, 1)));
% transposte deve ser usada no lugar do hermitiano
x2=sh2+n; % sinal desejado na entrada do equalizador + ruido
 
eq_out = filter(w, 1, x2);
dataSymbolsOut = qamdemod(eq_out,M); 
dataOutMatrix = de2bi(dataSymbolsOut,k);
dataOut = dataOutMatrix(:); 

[numErrors,ber] = biterr(dataIn,dataOut);

%ber
figure(3)
plot3(real(eq_out),imag(eq_out),1:rep,'r.'); % plota da modulaÃ§Ã£o 16 QAM
figure(4)
plot(1:5000,real(erro.*erro),'red');
% axis([-1 600 -1 1]);
xlabel('plot(e.*e)');

