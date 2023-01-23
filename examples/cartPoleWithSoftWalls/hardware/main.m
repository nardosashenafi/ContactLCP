binSize = 3;
binLayerSize = [5, 4, 3];
controlLayerSize = cell(binSize,1);

for i =1:binSize
    controlLayerSize{i} = [8, 4, 1];
end

binTotalParam = 69;
controlTotalParam = [89, 89,89];
param = rand(336, 1);

moe = MOE(binSize, binLayerSize, controlLayerSize, binTotalParam, controlTotalParam, param);
% x = [0.0, 3.14, 0.0, 0.2];
% 
% moe.control(x)