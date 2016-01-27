
function y = slicer_qpsk(u); 
 
% This function maps its complex input u 
% into the closest entry in a QPSK constellation.  

if real(u) >= 0   % Avoid using the sign function of 
   sign_a = 1;    % matlab since it returns 0 when its input
else              % is zero. This way y is always guaranteed
   sign_a = -1;   % to belong to the QPSK constellation
end 
 
if imag(u) >= 0 
   sign_b = 1; 
else 
   sign_b = -1; 
end 


y = sign_a + j*sign_b;
