clc;clear;
%%
%data_num=test_num;
load opt_data_attractor3_test1.mat
load min_rng_set_attractor3_test1.mat
numCopies = 1;  % 需要复制的次数
testdata = repmat(testdata, numCopies, 1);
testnum=size(testdata,1)/13;
%%
result = getfield ( opt_trials,'Fval');
 param= getfield ( opt_trials,'X');
 [sort_result,result_num]=sort(result);
 sort_param=param(result_num,:);
 opt_result=sort_param(1,:);
 sort_rng=min_rng_set(result_num);
 opt_rng=sort_rng(1);  
 rng(opt_rng);
 %%
eig_rho =opt_result(1);
W_in_a =opt_result(2);
a = opt_result(3);
reg = opt_result(4);
density =opt_result(5);
%%
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
 % 初始化图形
figure;
hold on;
for i=0:testnum-1
    
    Y1(:,1)=rand(12,1);
    x1=2*rand(resSize,1)-1;
    
    for t = 1:test_Len-1 
        x1 = (1-a)*x1 + a*tanh( Win*u + W*x1 );
        y = Wout*[1;x1;x1.^2;];
        % y = Wout*x1;
        Y1(:,t+1) = y;
        u(1:12) = y;  
        u(13)=testdata(13*(i+1),1+t);       
    end
    
    % 不同的 i 画不同的维度的相图
    if mod((i+1),3)==0
        save Loz.mat Y1 testdata
        % i+1 为 3 的倍数时绘制三维图 (维度 10, 11, 12)
        plot3(Y1(10,1000:testLen),Y1(11,1000:testLen),Y1(12,1000:testLen),'b-');
        plot3(testdata(13*i+10,1000:testLen),testdata(13*i+11,1000:testLen),testdata(13*i+12,1000:testLen),'r-');
    elseif mod((i+1),2)==0
        save 4D_ros.mat Y1 testdata
        % i+1 为 2 的倍数时绘制三维图 (维度 10, 11, 12)
        plot3(Y1(10,200:testLen),Y1(11,200:testLen),Y1(12,200:testLen),'b-');
        plot3(testdata(13*i+10,1000:testLen),testdata(13*i+11,1000:testLen),testdata(13*i+12,1000:testLen),'r-');
    else
        % 其他情况，绘制二维图 (维度 11, 12)
        save FHN.mat Y1 testdata
        plot(Y1(11,1000:testLen),Y1(12,1000:testLen),'b-');
        plot(testdata(13*i+11,1000:testLen),testdata(13*i+12,1000:testLen),'r-');
    end
    
end

rmse_dynamic=rmse_dynamic+mean(abs(Y1(11,300:testLen)-testdata(13*i+11,300:testLen)))+...
    mean(abs(Y1(12,300:testLen)-testdata(13*i+12,300:testLen)));
end
rmse_dynamic=rmse_dynamic/testnum;
rmse=rmse_dynamic;
if isnan(rmse) || rmse>10
    rmse=10;
end

