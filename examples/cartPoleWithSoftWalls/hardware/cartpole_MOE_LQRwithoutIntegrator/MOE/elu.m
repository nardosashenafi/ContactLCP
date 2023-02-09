function a = elu(x)
    a = zeros(length(x), 1);
    for i = 1:length(x)
         if x(i) > 0 
            a(i) = x(i); 
         else
            a(i) = 1.0 .* (exp(x(i)) - 1);
         end
    end
end