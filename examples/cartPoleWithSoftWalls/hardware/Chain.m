classdef Chain 
    properties 
        W1 = zeros(8, 5);
        W2 = zeros(4, 8);
        W3 = zeros(1, 4);
        b1 = zeros(5, 1);
        b2 = zeros(4, 1);
        b3 = zeros(8, 1);
    end
    
    methods
        
        function this = Chain(layerSize, vec)
            
            inputNum=5;
            layerNum = length(layerSize);
            layers = cat(2, inputNum, layerSize);
            layerRanges = zeros(layerNum, 4);
            % parse the weights and biases from the vector of parameters
            wbegins = 1;
            for i=2:length(layers)
               wends = wbegins + layers(i-1)*layers(i) - 1;
               bbegins = wends+1;
               bends = wends + layers(i);
               layerRanges(i-1,:) = [wbegins, wends, bbegins, bends];
               wbegins = bends + 1;
            end
            
            this.W1 = reshape(vec(layerRanges(1, 1):layerRanges(1, 2)), layers(2), layers(1));
            this.W2 = reshape(vec(layerRanges(2, 1):layerRanges(2, 2)), layers(3), layers(2));
            this.W3 = reshape(vec(layerRanges(3, 1):layerRanges(3, 2)), layers(4), layers(3));

            this.b1 = reshape(vec(layerRanges(1, 3):layerRanges(1, 4)), layers(2), 1);
            this.b2 = reshape(vec(layerRanges(2, 3):layerRanges(2, 4)), layers(3), 1);
            this.b3 = reshape(vec(layerRanges(3, 3):layerRanges(3, 4)), layers(4), 1);
        end
        

        function this = set.W1(this, value)
           this.W1 = value; 
        end
        function this = set.W2(this, value)
           this.W2 = value; 
        end

        function this = set.W3(this, value)
           this.W3 = value; 
        end

        function this = set.b1(this, value)
           this.b1 = value; 
        end

        function this = set.b2(this, value)
           this.b2 = value; 
        end

        function this = set.b3(this, value)
           this.b3 = value; 
        end
        
        function a = get.W1(this)
           a = this.W1; 
        end
        function a = get.W2(this)
           a = this.W2; 
        end

        function a = get.W3(this)
           a = this.W3; 
        end

        function a = get.b1(this)
           a = this.b1; 
        end

        function a = get.b2(this)
           a = this.b2; 
        end

        function a = get.b3(this)
           a = this.b3; 
        end
        
        function a = elu(this, x)
             if x >= 0 
                a = x; 
             else
                a = 1.0 .* (exp(x) - 1);
             end
        end
        
        function y = forward(this, x)
            z0 = x;
            z1 = this.elu(this.W1*z0 + this.b1);
            z2 = this.elu(this.W2*z1 + this.b2);
            z3 = this.elu(this.W3*z2 + this.b3);
            y = z3;
        end
        
        function y = forwardSoftmax(this, x)
            y = softmax(this.forward(x));
        end
    end
end