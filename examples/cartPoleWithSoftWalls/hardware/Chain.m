classdef Chain
    properties
        inputNum = 5;
        layerNum = 3;
        layerSize = [5, 4, 1];  
        W = cell(3, 1);
        b = cell(3, 1);
    end
    
    methods
        
        function this = Chain(layerSize, vec)
            this.layerSize = layerSize;
            this.layerNum = length(this.layerSize);
            layers = cat(2, this.inputNum, this.layerSize);
            wbegins = 1;
            for i=2:length(layers)
               wends = wbegins + layers(i-1)*layers(i) - 1;
               bsize = wends+1 : wends + layers(i);
               this.W{i-1} = reshape(vec(wbegins:wends),[layers(i), layers(i-1)]);
               this.b{i-1} = vec(bsize);
               wbegins = bsize(end) + 1;
            end
            
        end
        
        function a = elu(this, x)
             if x >= 0 
                a = x; 
             else
                a = 1*(exp(x) - 1);
             end
        end
        
        function y = forward(this, x)
            z = x;
            for i=1:this.layerNum
                z = this.elu(this.W{i}*z + this.b{i});
            end
            y = z;
        end
        
        function y = forwardSoftmax(this, x)
            y = softmax(this.forward(x));
        end
    end
end