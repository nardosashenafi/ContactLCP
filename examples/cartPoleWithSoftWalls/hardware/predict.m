inputNum = 5;
layerNum = 3;
biasConnect = [1;1;1];
inputConnect = [1 1 1 1 1 ;0 0 0 0 0 ;0 0 0 0 0];
layerConnect = [0 0 0; 1 0 0; 0 1 0];
outputConnect = [0 0 1];
net = network(inputNum, layerNum, biasConnect, inputConnect,layerConnect,outputConnect);
%creates one input layer and 2 additional layers, where the second layer is the
%output layer.
%input -> IW -> layer{1} -> LW{2,1} -> layer{2}

%set parameters size 
net.layers{1}.dimensions = 8;
net.layers{2}.dimensions = 4;
net.layers{3}.dimensions = 1;

%initialize parameters
net.inputs{1}.size = 1;
net.inputs{2}.size = 1;
net.inputs{3}.size = 1;
net.inputs{4}.size = 1;
net.inputs{5}.size = 1;
net.IW{1,1} = rand(8, 1);
net.IW{1,2} = rand(8, 1);
net.IW{1,3} = rand(8, 1);
net.IW{1,4} = rand(8, 1);
net.IW{1,5} = rand(8, 1);

net.LW{2,1} = rand(4, 8);   
net.LW{3,2} = rand(1, 4);
net.b{1} = rand(8, 1);     %2
net.b{2} = rand(4, 1);     %2
%net.b{3} = rand(2, 1);     %1
net.b{3} = rand(1);        %1
input = rand(inputNum, 1)

net(input)