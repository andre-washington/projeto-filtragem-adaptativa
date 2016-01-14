function [Rest] = estimar_matriz_autocorrelacao()

% Número de amostras para estimação
amostras = 500000;

% Canal
H = [1; 1.6];

% Dados originais
data = randn(amostras,1);

% Saída do canal
x = conv(data,H);
R = zeros(2,2);

for i = 2:length(x)
    v = [x(i); x(i - 1)];
    Raux = v * v';
    R = R + Raux;
end

Rest = R/(length(x)-1);