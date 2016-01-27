
function y = slicer16(u); 
 
% This function maps its complex input u 
% into the closest entry in a 16QAM constellation.  


if abs(real(u)) <= 2
   a = 1; 
else 
   a = 3;
end 
 
if abs(imag(u)) <= 2
   b = 1; 
else 
   b = 3;
end 


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


y = a*sign_a + j*b*sign_b;
