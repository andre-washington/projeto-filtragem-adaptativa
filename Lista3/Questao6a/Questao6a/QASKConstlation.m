function [xt, yt] = QASKConstlation(K) 
% Output the constellation in-phase and quadrature components. 
 
xx = constlay(K, 1); 
 
[leny, lenx] = size(xx); 
[xt, yt]= meshgrid([1-lenx : 2 : lenx-1], [leny-1 : -2 : 1-leny]'); 
 
xx = xx(:); 
xt = xt(:); 
yt = yt(:); 
 
tmp = isnan(xx); 
if ~isempty(tmp) 
  xx(tmp) = []; 
  xt(tmp) = []; 
  yt(tmp) = []; 
end; 
 
[xx, tmp] = sort(xx); 
xt = xt(tmp); 
yt = yt(tmp); 
end