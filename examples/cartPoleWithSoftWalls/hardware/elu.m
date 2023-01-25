function a = elu(x)
     if x >= 0 
        a = x; 
     else
        a = 1.0 .* (exp(x) - 1);
     end
end