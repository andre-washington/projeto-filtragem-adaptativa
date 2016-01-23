% Lista 3 - Questão 04
% Aluno: Davyd Bandeira de Melo

% Limpa o workspace, o console e fecha todas as janelas
clear, clc, close all;

% Número de amostras
amostras = 500;

% Coeficientes da função de transferência do canal
H = [1; 1.6];

% Dados gaussianos
dados = randn(amostras,1);

% Aplicação dos dados de entrada a função de transferência do sistema
% através de convolução
x = conv(dados,H);

% Inicialização dos coeficientes dos filtros para cada algoritmo
wgd        = [-2;-1]; % Gradiente deterministico
wn         = [-2;-1]; % Algoritmo de Newton
wlms       = [-2;-1]; % LMS
wlmsn      = [-2;-1]; % LMS Normalizado

% Valores anteriores dos pesos
wgd_ant    = wgd;
wn_ant     = wn;
wlms_ant   = wlms;
wlmsn_ant  = wlmsn;

% Realiza a estimação da matriz de correlação média
R = estimar_matriz_autocorrelacao();

% Realiza a estimação do vetor de correlação cruzada médio
p = estimar_correlacao_cruzada();

% Calcula a superfície de erro
m = -2:0.1:2;
n = m;
J = zeros(length(m),length(n));
for i=1:length(m)
    for j=1:length(n)
        aux = [m(i);n(j)];
        J(i,j) = 1 - 2*aux'*p + aux'*R*aux;
    end
end
J = J';

% Exibição dos gráficos
%main_plots = figure;
subplot(2,2,1);
contour(m,n,J,8);
title('Gradiente Determin�stico');
subplot(2,2,2);
contour(m,n,J,8);
title('Algorítmo de Newton');
subplot(2,2,3);
contour(m,n,J,8);
title('LMS');
subplot(2,2,4);
contour(m,n,J,8);
title('LMS Normalizado');

%EMQ
emq_GD      = zeros(amostras-1,1);
emq_N       = zeros(amostras-1,1);
emq_LMS     = zeros(amostras-1,1);
emq_LMSN    = zeros(amostras-1,1);

for i=2:amostras
    %Vetor acumulador
    x_v = [x(i); x(i-1)];
    
    %Sa�da desejada
    d = dados(i-1);
    
    %C�lculo do erro
    egd    = d - wgd'   * x_v;
    en     = d - wn'    * x_v;
    elms   = d - wlms'  * x_v;
    elmsn  = d - wlmsn' * x_v;
    
    wgd = wgd - 2 * 0.01 * (R * wgd - p);                          % Gradiente determinístico
    wn = wn - 0.1 * inv(R) * (-2 * p + 2 * R * wn);                % Algoritmo de Newton
    wlms = wlms + 2 * 0.01 * elms * x_v;                           % LMS
    wlmsn = wlmsn + 0.02 * (1/(x_v' * x_v + 0.01)) * elmsn .* x_v; % LMSN
    
    % Evolução do filtro
    
    %Gradiente Determin�stico
    subplot(2,2,1); hold on;
    plot([wgd_ant(1) wgd(1)], [wgd_ant(2) wgd(2)]);
    wgd_ant = wgd;
    hold off;
    
    %M�todo de Newton
    subplot(2,2,2); hold on;
    plot([wn_ant(1) wn(1)], [wn_ant(2) wn(2)]);
    wn_ant = wn;
    hold off;
    
    %LMS
    subplot(2,2,3); hold on;
    plot([wlms_ant(1) wlms(1)], [wlms_ant(2) wlms(2)]);
    wlms_ant = wlms;
    hold off;
    
    %LMSN
    subplot(2,2,4); hold on;
    plot([wlmsn_ant(1) wlmsn(1)], [wlmsn_ant(2) wlmsn(2)]);
    wlmsn_ant = wlmsn;
    hold off;
    
    % Erro
    emqgd(i-1)   = egd   * egd;
    emqn(i-1)    = en    * en;
    emqlms(i-1)  = elms  * elms;
    emqlmsn(i-1) = elmsn * elmsn;
end

w0 = inv(R)*p;

subplot(2,2,1); hold on;
plot(w0(1),w0(2),'r*');
hold off;

subplot(2,2,2); hold on;
plot(w0(1),w0(2),'r*');
hold off;

subplot(2,2,3); hold on;
plot(w0(1),w0(2),'r*');
hold off;

subplot(2,2,4); hold on;
plot(w0(1),w0(2),'r*');
hold off;

figure;
subplot(2,2,1);
plot(emqgd);
title('EMQ do Gradiente Determinístico'); xlabel('Passo'); ylabel('EMQ');

subplot(2,2,2);
plot(emqn);
title('EMQ do Algorítmo de Newton'); xlabel('Passo'); ylabel('EMQ');

subplot(2,2,3);
plot(emqlms);
title('EMQ do LMS'); xlabel('Passo'); ylabel('EMQ');
subplot(2,2,4);

plot(emqlmsn);
title('EMQ do LMS normalizado'); xlabel('Passo'); ylabel('EMQ');