 
function y = slicer64(u);  
  
% This function maps its complex input u  
% into the closest entry in a 64QAM constellation.   


if abs(real(u)) <= 2
   a = 1;  
elseif abs(real(u)) <= 4
   a = 3; 
elseif abs(real(u)) <= 6
   a = 5; 
else  
   a = 7;  
end  
  
if abs(imag(u)) <= 2
   b = 1;  
elseif abs(imag(u)) <= 4
   b = 3; 
elseif abs(imag(u)) <= 6
   b = 5; 
else  
   b = 7;  
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
