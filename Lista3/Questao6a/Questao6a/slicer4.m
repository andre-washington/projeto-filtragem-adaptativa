
function y = slicer4(u); 
 
% This function maps its complex input u 
% into the closest entry in a 4QAM constellation.  


if real(u) >= 0  
   sign_a = 1;   
else             
   sign_a = -1;  
end 
 
if imag(u) >= 0 
   sign_b = 1; 
else 
   sign_b = -1; 
end 


y = sign_a + j*sign_b;
