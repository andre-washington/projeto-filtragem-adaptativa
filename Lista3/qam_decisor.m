function [ output ] = qam_decisor( input, const_size )
switch(const_size)
    case 4
        if(abs(real(input)) == 0)
               switch(randi(2))
                    case 1
                        a = 1;
                    case 2
                        a = -1;
               end
        else
            a = sign(real(input));
        end
            
         if(abs(imag(input)) == 0)
               switch(randi(2))
                    case 1
                        b = 1;
                    case 2
                        b = -1;
               end
        else
            b = sign(imag(input));
        end    
    
    case 16
        if(abs(real(input)) == 0 && abs(real(input)) == 2)
              switch(randi(2))
                  case 1
                       a = real(input) + 1;
                  case 2
                       a = real(input) - 1;
              end
                
        else
               if (abs(real(input)) < 2)
                   a = 1*sign(real(input));
               else
                   a = 3*sign(real(input));
               end
        end
            %rounding imaginary part

        if(abs(imag(input)) == 0 && abs(imag(input)) == 2)
              switch(randi(2))
                  case 1
                       b = imag(input) + 1;
                  case 2
                       b = imag(input) - 1;
              end
                
        else
               if (abs(imag(input)) < 2)
                   b = 1*sign(imag(input));
               else
                   b = 3*sign(imag(input));
               end
        end
        
    case 64
        if(abs(real(input)) == 0 && abs(real(input)) == 2 ...
                && abs(real(input)) == 4 && abs(real(input)) == 6)
              switch(randi(2))
                  case 1
                       a = real(input) + 1;
                  case 2
                       a = real(input) - 1;
              end
              
        else
            if (abs(real(input)) < 2)
                a = 1*sign(real(input));
            else
                if (abs(real(input)) < 4)
                    a = 3*sign(real(input));
                else
                    if (abs(real(input)) < 6)
                        a = 5*sign(real(input));
                    else
                        a = 7*sign(real(input)); 
                    end
                end
            end
        end
            %rounding imaginary part

            if(abs(imag(input)) == 0 && abs(imag(input)) == 2 ...
                && abs(imag(input)) == 4 && abs(imag(input)) == 6)
              switch(randi(2))
                  case 1
                       b = imag(input) + 1;
                  case 2
                       b = imag(input) - 1;
              end
            else
                if (abs(imag(input)) < 2)
                        b = 1*sign(imag(input));
                else
                    if (abs(imag(input)) < 4)
                        b = 3*sign(imag(input));
                    else
                         if (abs(imag(input)) < 6)
                             b = 5*sign(imag(input));
                         else
                             b = 7*sign(imag(input)); 
                         end
                     end
                end
            end
        
    case 256
        if(abs(real(input)) == 0 && abs(real(input)) == 2 ...
                && abs(real(input)) == 4 && abs(real(input)) == 6 ...
                && abs(real(input)) == 8 && abs(real(input)) == 10 ...
                && abs(real(input))== 12 && abs(real(input)) == 14)
             switch(randi(2))
                  case 1
                       a = real(input) + 1;
                  case 2
                       a = real(input) - 1;
             end
        else       
            if (abs(real(input)) < 2)
                a = 1*sign(real(input));
            else
                if (abs(real(input)) < 4)
                    a = 3*sign(real(input));
                else
                    if (abs(real(input)) < 6)
                        a = 5*sign(real(input));
                    else
                        if (abs(real(input)) < 8)
                            a = 7*sign(real(input));
                        else
                            if (abs(real(input)) < 10)
                             a = 9*sign(real(input));
                            else
                                if (abs(real(input)) < 12)
                                     a = 11*sign(real(input));
                                else
                                    if (abs(real(input)) < 14)
                                        a = 13*sign(real(input));
                                    else
                                        a = 15*sign(real(input));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        if(abs(imag(input)) == 0 && abs(imag(input)) == 2 ...
                && abs(imag(input)) == 4 && abs(imag(input)) == 6 ...
                && abs(imag(input)) == 8 && abs(imag(input)) == 10 ...
                && abs(imag(input))== 12 && abs(imag(input)) == 14)
             switch(randi(2))
                  case 1
                       b = imag(input) + 1;
                  case 2
                       b = imag(input) - 1;
             end
        else       
            if (abs(imag(input)) < 2)
                b = 1*sign(imag(input));
            else
                if (abs(imag(input)) < 4)
                    b = 3*sign(imag(input));
                else
                    if (abs(imag(input)) < 6)
                        b = 5*sign(imag(input));
                    else
                        if (abs(imag(input)) < 8)
                            b = 7*sign(imag(input));
                        else
                            if (abs(imag(input)) < 10)
                             b = 9*sign(imag(input));
                            else
                                if (abs(imag(input)) < 12)
                                     b = 11*sign(imag(input));
                                else
                                    if (abs(imag(input)) < 14)
                                        b = 13*sign(imag(input));
                                    else
                                        b = 15*sign(imag(input));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
end
    output = complex(a,b);
end
