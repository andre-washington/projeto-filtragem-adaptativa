function [ out ] = decisor( xw, limit )
%DECISOR retorna o valor de xd arredondado, com o limite definido por limit
    
    if(real(xw) >0)
        xw= ceil(real(xw)) + imag(xw)*1i;
    else 
        xw= floor(real(xw)) + imag(xw)*1i;
    end
     if(imag(xw) >0)
        xw= real(xw) + ceil(imag(xw))*1i;
    else 
        xw= real(xw) + floor(imag(xw))*1i;
    end
    
    xd = xw;
    
    if(real(xd)>limit)
        xd=limit+imag(xd)*1i;
    end
    if(real(xd)<-limit)
        xd= -limit + imag(xd)*1i;
    end
    if(imag(xd)>limit)
        xd= real(xd) + limit*i;
    end
    
    if(imag(xd)<-limit)
       xd= real(xd) - limit*i;
    end
    
    out = xd;
end

