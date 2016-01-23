% Lista 4 - Questão 04

% Limpa o workspace, o console e fecha todas as janelas
clear, clc, close all;

% Número de amostras
amostras = 100;

% Coeficientes da função de transferência do canal
H = [1; 1.6];

% Dados gaussianos
dados = randn(amostras,1);
%dados = sin(1:amostras)';
% Aplicação dos dados de entrada a função de transferência do sistema
% através de convolução
x = conv(dados,H);


% Realiza a estimação da matriz de correlação média
%R = estimar_matriz_autocorrelacao();
R =zeros(2,2);
xr=conv(dados,H);
for i=2:length(xr)
    v=[xr(i);xr(i-1)];
    R = R + v*v';
end
R = R/(amostras);
% Realiza a estimação do vetor de correlação cruzada médio
%p = estimar_correlacao_cruzada();
xp=conv(dados,H);
p=zeros(2,1);
for i=1:length(xp)-1
    v=[xp(i);xp(i+1)];
    p = p + v * dados(i);
end

% Calcula a superfície de erro
m = -10:0.1:10;
n = m;
J = zeros(length(m),length(n));
for i=1:length(m)
    for j=1:length(n)
        aux = [m(i);n(j)];
        J(i,j) = 1 - 2*aux'*p + aux'*R*aux;
    end
end
J = J';

% RLS ALTERNATIVO
num_ite = amostras;
num_taps = 2;
input = x(1:num_ite+1);
lambda = 0.3;
delta = 100;
w = [-10 10]';
sd = delta*eye(num_taps);
pd = zeros(num_taps, 1);

xRLS = zeros(num_taps, 1);
err_vec = zeros(num_ite, 1);
y = zeros(num_ite, 1);
for k = 1:num_ite 
    wplot(:,k)=w;
   d = input(k+1);
   if (k < num_taps)
       xRLS = [input(k); xRLS(1:end-1)];
   else
       xRLS = input(k:-1:k-num_taps+1);
   end
   e = d - xRLS'*w;
   phi = sd*xRLS;
   
   sd = 1/lambda*(sd - (phi*phi')/(lambda + phi'*xRLS));
   
   w = w + e*sd*xRLS;
   
   y(k) = w'*xRLS;
   err_vec(k) = d - y(k);
end
figure(1)
plot(input(2:end), 'r');
hold on;
plot(y, 'b--');
hold off;
figure(2)
contour(m,n,J,8);
hold on
plot(wplot(1,1:end),wplot(2,1:end));
title('RLS Alternativo');
figure(3)
plot(err_vec.^2);
title(['Num de Coef = 2  \lambda = ' num2str(lambda)]);
