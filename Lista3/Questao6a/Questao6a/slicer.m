
function y = slicer(input); 
 
 % This function maps its complex input signal
 % into the closest entry in the constellation.  
 
u = input(1) + j*input(2);
const_type = input(3);
  
 switch const_type
 case 16
    y = slicer16(u);
 case 64
    y = slicer64(u);
 case 256
    y = slicer256(u);
 end
 