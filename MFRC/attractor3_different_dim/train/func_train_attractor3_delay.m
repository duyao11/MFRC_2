function [rmse,rmse_dynamic]  = func_train_attractor3_delay(hyperpara_set,train_data,rng_num,testdata,data_num,data_len)
%%
rng(rng_num);
%data_num=test_num;
numCopies = 1;  % 需要复制的次数
testdata = repmat(testdata, numCopies, 1);
testnum=size(testdata,1)/13;
%%
eig_rho =hyperpara_set(1);
W_in_a = hyperpara_set(2);
a = hyperpara_set(3);
reg = hyperpara_set(4);
density =hyperpara_set(5);
resSize =100; % size of the reservoir nodes;  
initLen = 100;
TrainLen=sum(data_len)-1;
test_Len =2000;
testLen=2000;
inSize = 13; 
outSize = 12;
nonliner_num=2;
%%
  indata=train_data;
  X = zeros(nonliner_num*resSize+1,TrainLen);
  %X = zeros(resSize,TrainLen);
  Yt = indata(1:outSize,2:TrainLen+1);% run the reservoir with the data and collect X
%%
Win = (2.0*rand(resSize,inSize)-1.0)*W_in_a;
WW = zeros(resSize,resSize);
for i=1:resSize
    for j=i:resSize
            if (rand()<density)
             WW(i,j)=(2.0*rand()-1.0);
             WW(j,i)=WW(i,j);
            end
    end
end
rhoW = eigs(WW,1);
W = WW .* (eig_rho /rhoW); 
x = zeros(resSize,1);
for t = 1:TrainLen
    u = indata(:,t);
    x = (1-a)*x + a*tanh( Win*u + W*x );
    X(:,t) = [1;x;x.^2;];
    %X(:,t) = x;
end

Data_len=[0,data_len];
for i=0:data_num-1
   trainLen=sum(Data_len(1:i+1));
   X(:,trainLen+1-i*initLen:trainLen+initLen-i*initLen)=[];
   Yt(:,trainLen+1-i*initLen:trainLen+initLen-i*initLen)=[];
end
rank=randperm( size(X,2) );  
X=X(:, rank); 
Yt=Yt(:, rank); 
X_T = X';
Wout = Yt*X_T / (X*X_T + reg*eye(nonliner_num*resSize+1));
%%
Y1= zeros(outSize,testLen);
rmse_dynamic=0;  
rng(4);
for i=0:testnum-1
    
    
    Y1(:,1)=rand(12,1);
    x1=2*rand(resSize,1)-1;
for t = 1:test_Len-1 
   
    x1 = (1-a)*x1 + a*tanh( Win*u + W*x1 );
    y = Wout*[1;x1;x1.^2;];
    %y = Wout*x1;
    Y1(:,t+1) = y;
    u(1:12) = y;  
   
    u(13)=testdata(13*(i+1),1+t);
   
        
end
rmse_dynamic=rmse_dynamic+mean(abs(Y1(11,300:testLen)-testdata(13*i+11,300:testLen)))+...
    mean(abs(Y1(12,300:testLen)-testdata(13*i+12,300:testLen)));
end
rmse_dynamic=rmse_dynamic/testnum;
rmse=rmse_dynamic;
if isnan(rmse) || rmse>10
    rmse=10;
end
