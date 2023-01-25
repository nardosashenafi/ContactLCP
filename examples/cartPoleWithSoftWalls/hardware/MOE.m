classdef MOE < matlab.System
    properties

    end
    
    properties(DiscreteState)

    end
    
    properties(Access = private)
        
    end
    
    methods
%         %constructors
        function this = MOE(bin, controlArray1, controlArray2, controlArray3)
            %%Neural network settings and parameters
%             binSize = 3;
%             binLayerSize = [5, 4, 3];
%             controlLayerSize = [[8, 4, 1]; [8, 4, 1]; [8, 4, 1]];
% 
%             binTotalParam = 69;
%             controlTotalParam = [89, 89,89];
%             param = rand(336, 1);
% 
%             %parse parameters
%             nnParamRanges = zeros(binSize, 2);
%             paramBegin = 1;
%             paramEnd = paramBegin + binTotalParam - 1;
%             psi = param(paramBegin:paramEnd);
%             paramBegin = paramEnd+1;
%             for i = 1:binSize
%                 paramEnd = paramBegin + controlTotalParam(i) - 1;
%                 nnParamRanges(i,:) = [paramBegin, paramEnd];
%                 paramBegin = paramEnd+1;
%             end
%             
%             %construct gating network
%             this.bin = Chain(binLayerSize, psi);
%             %construct controllers and populate parameters
%             this.controlArray1 = Chain(controlLayerSize(1,:), param(nnParamRanges(1, 1):nnParamRanges(1, 2)));
%             this.controlArray2 = Chain(controlLayerSize(2,:), param(nnParamRanges(2, 1):nnParamRanges(2, 2)));
%             this.controlArray3 = Chain(controlLayerSize(3,:), param(nnParamRanges(3, 1):nnParamRanges(3, 2)));
%             
        end
        
        function u = control(this, x)
            z = this.inputLayer(x);
            pk = this.bin.forwardSoftmax(z);
            [pkmax, k] = max(pk);
            fname = strcat('this.controlArray', int2str(eval('k')), '.forward(z)');
            u = eval(fname);
        end
        
        function z = inputLayer(this, x)
             z = [x(1); cos(x(2)); sin(x(2)); x(3); x(4)];
        end
        
       
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods(Access = protected)
        function setupImpl(obj)
            
        end

        function y = stepImpl(obj, x1, theta, v, omega)
            x = cat(2, x1(1), theta(1), v(1), omega(1));
%             y = obj.control(x);
            z = obj.inputLayer(x);
            W = obj.bin.W1;
            y = sum(W*z);
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
      
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end

