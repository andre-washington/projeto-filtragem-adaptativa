 
function y = slicer256(u);  
  
% This function maps its complex input u  
% into the closest entry in a 256QAM constellation.   


if abs(real(u)) <= 2
   a = 1;  
elseif abs(real(u)) <= 4 
   a = 3; 
elseif abs(real(u)) <= 6
   a = 5; 
elseif abs(real(u)) <= 8
   a = 7; 
elseif abs(real(u)) <= 10
   a = 9;   
elseif abs(real(u)) <= 12
   a = 11;   
elseif abs(real(u)) <= 14
   a = 13;    
else  
   a = 15;  
end  
  
if abs(imag(u)) <= 2  
   b = 1;  
elseif abs(imag(u)) <= 4
   b = 3; 
elseif abs(imag(u)) <= 6
   b = 5; 
elseif abs(imag(u)) <= 8
   b = 7; 
elseif abs(imag(u)) <= 10
   b = 9; 
elseif abs(imag(u)) <= 12
   b = 11; 
elseif abs(imag(u)) <= 14
   b = 13; 
else  
   b = 15;
end  
 
 
if real(u) >= 0   % Avoid using the sign function of  
   sign_a = 1;    % matlab since it returns 0 when its input 
else              % is zero.
   sign_a = -1;  
end  
  
if imag(u) >= 0  
   sign_b = 1;  
else  
   sign_b = -1;  
end  
 
 
y = a*sign_a + j*b*sign_b; 
