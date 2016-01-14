% Lista 4 - Questão 01
% Aluno: Davyd Bandeira de Melo

% Limpa o workspace, o console e fecha todas as janelas
clear, clc, close all;

% Inicialização dos parâmetros do algoritmo
amostras = 100;
w = [1;0;0];
want = w;
delta = 100;
Sd = delta*eye(3);
l = 0.98;
sd = 0.5;

% Geração dos dados
x = cos((pi/3)*(0:amostras-1));
xant = [0;0;0];

% Estimativa das estatísticas dos dados
R = estimar_matriz_autocorrelacao();
p = estimar_correlacao_cruzada();

% Geração da superfície de erro
m = -2:0.1:2;
n = m;
J = zeros(length(m),length(n));
for i=1:length(m)
    for j=1:length(n)
        aux = [1; m(i); n(j)];
        J(i,j) = sd - 2 * aux' * p + aux' * R * aux;
    end
end
J = J';

figure; subplot(2,2,3); contour(m,n,J,8); title('Curvas de Contorno');

mse = zeros(amostras-3,1);

for i = 3:amostras-1
    
    x_v = x(i:-1:i-2)';
    d = x(i+1);
    
    error = d - x_v' * w; % Priori
    
    psi = Sd * x_v;
    Sd = (1/l)*(Sd - (psi * psi') / (l + psi' * x_v));
    w = w + error * Sd *x_v;
    w(1) = 1;
    
    error_pos = d - x_v'*w; % Posteriori
    
    subplot(2,2,1); hold on; plot([i-3 i-2], [want' * xant w' * x_v]); hold off;
    subplot(2,2,3); hold on; plot([want(2) w(2)], [want(3) w(3)]); hold off;
    
    want = w; xant = x_v;
    
    % MSE
    mse(i-2) = error_pos*error_pos;
end

% Exibição dos gráficos
subplot(2,2,1); title('Predição');
xlabel('Passo');

subplot(2,2,3);
title('Pesos do filtro');
xlabel('Passo');

subplot(2,2,2); plot(3:(amostras-1), cos((pi/3)*(4:amostras)));

title('Sinal original');
xlabel('Passo');

subplot(2,2,4);
plot(mse);
title('MSE (a posteriori) do RLS');
xlabel('Passo');