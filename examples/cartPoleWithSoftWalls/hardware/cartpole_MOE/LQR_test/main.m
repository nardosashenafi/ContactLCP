% lqr_gains = [-0.707107;  25.8856; -3.51898; 4.49824];
lqr_gains = 1.0*[-0.7;  25.0; -4.0; 5.0];
Fmax = 4.5;
binSize = 3;
binLayerSize = [4, 3, 3];
controlLayerSize = [[10, 4, 1]; [10, 4, 1]; [10, 4, 1]];

binTotalParam = 51;
controlTotalParam = [109, 109, 109];
savedWeights = readtable("../saved_weights/hardwareParam_friction12_mediumpumping.csv");
param = savedWeights.param; 

%parse parameters for each neural network
nnParamRanges = zeros(binSize, 2);
paramBegin = 1;
paramEnd = paramBegin + binTotalParam - 1;
psi = param(paramBegin:paramEnd); 
paramBegin = paramEnd+1;
for i = 1:binSize
    paramEnd = paramBegin + controlTotalParam(i) - 1;
    nnParamRanges(i,:) = [paramBegin, paramEnd];
    paramBegin = paramEnd+1;
end

%%parse bin parameters
[binLayerRanges, binLayers] = parseParams(binLayerSize);
binW1 = reshape(psi(binLayerRanges(1, 1):binLayerRanges(1, 2)), binLayers(2), binLayers(1));
binW2 = reshape(psi(binLayerRanges(2, 1):binLayerRanges(2, 2)), binLayers(3), binLayers(2));
binW3 = reshape(psi(binLayerRanges(3, 1):binLayerRanges(3, 2)), binLayers(4), binLayers(3));

binb1 = reshape(psi(binLayerRanges(1, 3):binLayerRanges(1, 4)), binLayers(2), 1);
binb2 = reshape(psi(binLayerRanges(2, 3):binLayerRanges(2, 4)), binLayers(3), 1);
binb3 = reshape(psi(binLayerRanges(3, 3):binLayerRanges(3, 4)), binLayers(4), 1);
  
%%parse controlArray1 parameters
[control1LayerRanges, control1Layers] = parseParams(controlLayerSize(1,:));
theta1 = param(nnParamRanges(1,1):nnParamRanges(1,2));

control1W1 = reshape(theta1(control1LayerRanges(1, 1):control1LayerRanges(1, 2)), control1Layers(2), control1Layers(1));
control1W2 = reshape(theta1(control1LayerRanges(2, 1):control1LayerRanges(2, 2)), control1Layers(3), control1Layers(2));
control1W3 = reshape(theta1(control1LayerRanges(3, 1):control1LayerRanges(3, 2)), control1Layers(4), control1Layers(3));

control1b1 = reshape(theta1(control1LayerRanges(1, 3):control1LayerRanges(1, 4)), control1Layers(2), 1);
control1b2 = reshape(theta1(control1LayerRanges(2, 3):control1LayerRanges(2, 4)), control1Layers(3), 1);
control1b3 = reshape(theta1(control1LayerRanges(3, 3):control1LayerRanges(3, 4)), control1Layers(4), 1);

%%parse controlArray2 parameters
[control2LayerRanges, control2Layers] = parseParams(controlLayerSize(2,:));
theta2 = param(nnParamRanges(2,1):nnParamRanges(2,2));

control2W1 = reshape(theta2(control2LayerRanges(1, 1):control2LayerRanges(1, 2)), control2Layers(2), control2Layers(1));
control2W2 = reshape(theta2(control2LayerRanges(2, 1):control2LayerRanges(2, 2)), control2Layers(3), control2Layers(2));
control2W3 = reshape(theta2(control2LayerRanges(3, 1):control2LayerRanges(3, 2)), control2Layers(4), control2Layers(3));

control2b1 = reshape(theta2(control2LayerRanges(1, 3):control2LayerRanges(1, 4)), control2Layers(2), 1);
control2b2 = reshape(theta2(control2LayerRanges(2, 3):control2LayerRanges(2, 4)), control2Layers(3), 1);
control2b3 = reshape(theta2(control2LayerRanges(3, 3):control2LayerRanges(3, 4)), control2Layers(4), 1);

%%parse controlArray3 parameters
[control3LayerRanges, control3Layers] = parseParams(controlLayerSize(3,:));
theta3 = param(nnParamRanges(3,1):nnParamRanges(3,2));

control3W1 = reshape(theta3(control3LayerRanges(1, 1):control3LayerRanges(1, 2)), control3Layers(2), control3Layers(1));
control3W2 = reshape(theta3(control3LayerRanges(2, 1):control3LayerRanges(2, 2)), control3Layers(3), control3Layers(2));
control3W3 = reshape(theta3(control3LayerRanges(3, 1):control3LayerRanges(3, 2)), control3Layers(4), control3Layers(3));

control3b1 = reshape(theta3(control3LayerRanges(1, 3):control3LayerRanges(1, 4)), control3Layers(2), 1);
control3b2 = reshape(theta3(control3LayerRanges(2, 3):control3LayerRanges(2, 4)), control3Layers(3), 1);
control3b3 = reshape(theta3(control3LayerRanges(3, 3):control3LayerRanges(3, 4)), control3Layers(4), 1);