
clc;clear;
load opt_data_attractor4_range_0_beta_range_3.mat 
load min_rng_set_attractor4_range_0_beta_range_3.mat
%  load opt_data_attractor4_range4.mat 
%  load min_rng_set_attractor4_range4.mat
% load opt_data_attractor4_range_0.mat
% load min_rng_set_attractor4_range_0.mat

result = getfield ( opt_trials,'Fval');
 param= getfield ( opt_trials,'X');
 [sort_result,result_num]=sort(result);
 sort_param=param(result_num,:);
 opt_result=sort_param(1,:);
 sort_rng=min_rng_set(result_num);
 opt_rng=sort_rng(1);
 rng(opt_rng);
%rng(10);
drive_num=1;
nonliner_num=2;
numCopies =3;  % 需要复制的次数
testdata = repmat(testdata, numCopies, 1);
%Y0_set=0.8*(2*rand(6,2)-1);
%testdata=func_gen_test_data_chua_para(Y0_set,c3);

%Testdata=testdata;
%testdata=1.6*ones(3,4000);

%testdata=[testdata;Testdata(1:3,:)];
%testdata=[testdata;Testdata(28:30,:)];
testnum=size(testdata,1)/4;
%testnum=1;
%%
eig_rho =opt_result(1);
W_in_a =opt_result(2);
a = opt_result(3);
reg = opt_result(4);
density =opt_result(5);
resSize =100; % size of the reservoir nodes;  
initLen = 100;
TrainLen=sum(data_len)-1;
test_Len = 3000;
testLen=3000;
inSize = 4; 
outSize = 3;
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
% figure
% plot(Wout(1,2:101),'g.');
%%
rmse=0;
rmse_dynamic=0;    
color_f1= addcolorplus(103);
color_f2= addcolorplus(101);
color_f3= addcolorplus(13);
color_f4= addcolorplus(1);
color_f5=addcolorplus(118);
Y1= zeros(outSize,testLen);
figure
rng(4);
RC_loz=[];
RC_ros=[];
RC_lim=[];
RC_qua=[];
RC_out1=[];
RC_out2=[];
RC_out3=[];
j=0;
error=zeros(1,testnum);
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
    u = y;     
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
 RC_out1=[RC_out1,Y1(1,1200:end)];
 RC_out2=[RC_out2,Y1(2,1200:end)];
 RC_out3=[RC_out3,Y1(3,1200:   end)];
   error(i+1)=abs(cohens_d1)+abs(cohens_da)+abs(cohens_dm);
if mod(i+2,4)==0
%  figure
%  plot(Y1(2,1000:end),'Color',color_f2,'linewidth',1);
%  hold on;
%  set(gca,'FontName','Times New Roman','FontSize',18);
  figure
 plot3(Y1(1,1000:end),Y1(2,1000:end),Y1(3,1000:end),'Color',color_f2,'linewidth',1);
 hold on;
end
end

threlod=1;
error_num=0;
for i=0:numCopies-1
    if error(4*i+1)>threlod 
       error_num=error_num+1;
    end
end
accuracy=(numCopies-error_num)/numCopies;


figure
plot3(train_data(1,:),train_data(2,:),train_data(3,:),'linestyle','--');
set(gca,'FontName','Times New Roman','FontSize',18);
hold on;
plot3(RC_out1,RC_out2,RC_out3,'Color',color_f5,'linewidth',1);
set(gca,'FontName','Times New Roman','FontSize',18);


% 

% plot3(RC_out1(1:1799),RC_out2(1:1799),RC_out3(1:1799),'Color',color_f1,'linewidth',1);
% hold on;
% plot3(RC_out1(1802:3599),RC_out2(1802:3599),RC_out3(1802:3599),'Color',color_f4,'linewidth',1);
% hold on;
% plot3(RC_out1(3603:5399),RC_out2(3603:5399),RC_out3(3603:5399),'Color',color_f3,'linewidth',1);
% hold on;
% plot3(RC_out1(5413:7200),RC_out2(5413:7200),RC_out3(5413:7200),'Color',color_f2,'linewidth',1);
% set(gca,'FontName','Times New Roman','FontSize',18);

