function [rmse,rmse_dynamic]  = func_train_chua_no_para(hyperpara_set,train_data,rng_num,testdata,data_num,data_len)
%%
rng(rng_num);
drive_num=1;

%data_num=test_num;
numCopies = 3;  % 需要复制的次数
testdata = repmat(testdata, numCopies, 1);
testnum=size(testdata,1)/4;
%%
eig_rho =hyperpara_set(1);
W_in_a = hyperpara_set(2);
a = hyperpara_set(3);
reg = hyperpara_set(4);
density =hyperpara_set(5);
resSize =100; % size of the reservoir nodes;  
initLen = 100;
TrainLen=sum(data_len)-1;
test_Len = 3000;
testLen=3000;
inSize = 4; 
outSize = 3;
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
    
    
    Y1(:,1)=rand(3,1);
for t = 1:test_Len-1 
   if t<=drive_num
    x1=2*rand(resSize,1)-1;
    else
    x1 = (1-a)*x1 + a*tanh( Win*u + W*x1 );
   end
    y = Wout*[1;x1;x1.^2;];
    %y = Wout*x1;
    Y1(:,t+1) = y;
     if t<drive_num
    u(1,1)=testdata(4*(i+1)-3,1+t);
     u(2,1)=testdata(4*(i+1)-2,1+t);
     u(3,1)=testdata(4*(i+1)-1,1+t);
    else
    u = y;   
     end  
    u(4,1)=testdata(4*(i+1),1+t);
end
group1=testdata(4*i+1,1:testLen);
group2=Y1(1,:);
groupa=testdata(4*i+2,1:testLen);
groupb=Y1(2,:);
groupm=testdata(4*i+3,1:testLen);
groupn=Y1(3,:);

% 计算均值和标准差
mean1 = mean(group1);
mean2 = mean(group2);
meana = mean(groupa);
meanb = mean(groupb);
meanm = mean(groupm);
meann = mean(groupn);

std1 = std(group1);
std2 = std(group2);
stda = std(groupa);
stdb = std(groupb);
stdm = std(groupm);
stdn = std(groupn);


% 计算合并标准差
pooled_std1 = sqrt(((length(group1) - 1) * std1^2 + (length(group2) - 1) * std2^2) / (length(group1) + length(group2) - 2));

% 计算Cohen's d
cohens_d1 = (mean1 - mean2) / pooled_std1;
 % 计算合并标准差
pooled_stda = sqrt(((length(groupa) - 1) * stda^2 + (length(groupb) - 1) * stdb^2) / (length(groupa) + length(groupb) - 2));

% 计算Cohen's d
cohens_da = (meana - meanb) / pooled_stda;
% 计算合并标准差
pooled_stdm = sqrt(((length(groupm) - 1) * stdm^2 + (length(groupn) - 1) * stdn^2) / (length(groupm) + length(groupn) - 2));
% 计算Cohen's d
cohens_dm = (meanm - meann) / pooled_stdm;
rmse_dynamic=rmse_dynamic+abs(cohens_d1)+abs(cohens_da)+abs(cohens_dm);
end
rmse_dynamic=rmse_dynamic/testnum;
rmse=rmse_dynamic;
if isnan(rmse) || rmse>10
    rmse=10;
end
