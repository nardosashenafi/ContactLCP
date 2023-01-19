inputNum = 4;
layerNum = 3;
biasConnect = [1;1;1];
inputConnect = [1 1 1 1;0 0 0 0;0 0 0 0];
layerConnect = [0 0 0; 1 0 0; 0 1 0];
outputConnect = [0 0 1];
net = network(inputNum, layerNum, biasConnect, inputConnect,layerConnect,outputConnect);

%set parameters size 
net.layers{1}.dimensions = 8;
net.layers{2}.dimensions = 4;
net.layers{3}.dimensions = 2;

%initialize parameters
net.LW{2,1} = rand(2, 4)    %(2, 4)
net.LW{3, 2} = rand(1, 2)   %(1, 2)
net.biases{2} = rand(2)     %2
net.biases{3} = rand(1)     %1