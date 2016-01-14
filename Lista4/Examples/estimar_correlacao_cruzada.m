function [pxdest] = estimar_correlacao_cruzada()

% Número de amostras para estimação
amostras = 500000;

% Canal
H = [1; 1.6];

% Dados originais
data = randn(amostras,1);

% Saída do canal
x = conv(data,H);
pxd = zeros(2,1);

for i = 1:length(x)-1
    v = [x(i); x(i + 1)];
    pxd_aux = v * data(i);
    pxd = pxd + pxd_aux;
end

pxdest = pxd/(length(x)-2);
pxdest = pxdest(end:-1:1);