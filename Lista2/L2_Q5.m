% Filtragem Adaptativa: Filtro ótimo. 
% questão 5 lista 2
clc;
clear all;
dsigma = 24.4;
w0 = -20:0.05:20;
w1 = -20:0.05:20;
% [x,y] = meshgrid(w0,w1);
jw=zeros(length(w0),length(w1));
for i = 1:length(w0)
    for j=1:length(w1)
        jw(i,j) = w0(i)^2 + w1(j)^2 - 4*w0(i) - 9*w1(j) + dsigma;
    end
end
figure(1);
mesh (w0, w1,jw)
 figure(2);
zmin = floor(min(jw(:)));
zmax = ceil(max(jw(:)));
zinc = (zmax - zmin) / 40;
zlevs = zmin:zinc:zmax;

contour(w0,w1,jw,zlevs) 
axis([-25 25 -25 25])