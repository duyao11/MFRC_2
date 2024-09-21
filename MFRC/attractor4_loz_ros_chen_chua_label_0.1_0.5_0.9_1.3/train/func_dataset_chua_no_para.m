%clc;clear;
function [train_data,testdata,data_num,data_len]= func_dataset_chua_no_para()
data_num=4;
data_len=ones(1,data_num)*3000;
range=0.2;
range_lim=[-range,range,range];
range_qua=[range,-range,-range];
range_chen=[-range,range,range];
range_chua=[range,-range,-range];
% range_chen=[-range,-range,range];
% range_chua=[-range,range,-range];
mod_num=20;
transient=40002;
%% Lorenz
F=@(t,Y)[10*Y(2)-10*Y(1);...
         -Y(1).*Y(3)+28*Y(1)-Y(2);...
         Y(1).*Y(2)-8*Y(3)/3;...
         ];
         
opt = odeset( 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));     
y0=rand(3,1);
[t,Y] = ode45(F, [0:0.001:100], y0, opt);
data_lorenz=Y(transient:end,:);
max_loz=zeros(3,1);
min_loz=zeros(3,1);
for i=1:3
max_loz(i) = max(data_lorenz(:,i));
min_loz(i) = min(data_lorenz(:,i));
end
normalized_loz=data_lorenz';
% 归一化
data_loz=[];
for i=1:3
normalized_loz(i,:) = (data_lorenz(:,i)' - min_loz(i)) / (max_loz(i)-min_loz(i))/2+range;
xf1=normalized_loz(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num)==0);
data_loz=[data_loz;xf1];
end

%% Rossler
F=@(t,Y)[-5*Y(2)-5*Y(3);...
         5*Y(1)+5*Y(2)/2;...
         10+5*Y(3)*(Y(1)-4);...
         ];   
opt = odeset( 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));     
y0=rand(3,1);
[t,Y] = ode45(F, [0:0.001:100], y0, opt);
data_rossler=Y(transient:end,:);
max_rossler=zeros(3,1);
min_rossler=zeros(3,1);
for i=1:3
max_rossler(i) = max(data_rossler(:,i));
min_rossler(i) = min(data_rossler(:,i));
end
normalized_rossler=data_rossler';
% 归一化
data_ros=[];
for i=1:3
normalized_rossler(i,:) = (data_rossler(:,i)' - min_rossler(i)) /(max_rossler(i)- min_rossler(i))/2-range;
xf1=normalized_rossler(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num)==0);
data_ros=[data_ros;xf1];
end

%% limit

F=@(t,Y)[10*Y(1)*(2-Y(1).*Y(1)-Y(2).*Y(2))-10*Y(2);...
         10*Y(2).*(2-Y(1).*Y(1)-Y(2).*Y(2))+10*Y(1);...
         -10*Y(3);...
         ];   
opt = odeset( 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));     
y0=rand(3,1);
[t,Y] = ode45(F, [0:0.001:100], y0, opt);
data_limit=Y(transient:end,:);
max_limit=zeros(3,1);
min_limit=zeros(3,1);
for i=1:3
max_limit(i) = max(data_limit(:,i));
min_limit(i) = min(data_limit(:,i));
end
%min_limit(3)=1;
normalized_limit=data_limit';
% 归一化
data_lim=[];
for i=1:3
normalized_limit(i,:) = (data_limit(:,i)' - min_limit(i)) / (max_limit(i)-min_limit(i))/2+range_lim(i);
xf1=normalized_limit(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num)==0);
data_lim=[data_lim;xf1];
end
%% quasi 
F=@(t,Y)[2*pi*Y(2);...
        2*pi*Y(2).*(Y(1).*Y(1)-0.5*Y(1).^4)-4*pi^2*Y(1);...
         0.9-Y(1).*Y(1);...
         ];   
opt = odeset( 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));     
y0=rand(3,1);
[t,Y] = ode45(F, [0:0.001:100], y0, opt);
data_quasi=Y(transient:end,:);
max_quasi=zeros(3,1);
min_quasi=zeros(3,1);
for i=1:3
max_quasi(i) = max(data_quasi(:,i));
min_quasi(i) = min(data_quasi(:,i));
end
normalized_quasi=data_quasi';
% 归一化
data_qua=[];
for i=1:3
normalized_quasi(i,:) = (data_quasi(:,i)' - min_quasi(i)) / (max_quasi(i)-min_quasi(i))/2+range_qua(i);
xf1=normalized_quasi(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num)==0);
data_qua=[data_qua;xf1];
end
%%
%% Chen
a=40;b=3;c=28;
F=@(t,Y)[a*(Y(2)-Y(1));...
         (c-a)*Y(1)-Y(1)*Y(3)+c*Y(2);...
         Y(1)*Y(2)-b*Y(3);...
         ];        
opt = odeset( 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));     
y0=rand(3,1);
[t,Y] = ode45(F, [0:0.001:100], y0, opt);
data_chen=Y(transient:end,:);
max_chen=zeros(3,1);
min_chen=zeros(3,1);
for i=1:3
max_chen(i) = max(data_chen(:,i));
min_chen(i) = min(data_chen(:,i));
end
normalized_chen=data_chen';
% 归一化
data_Chen=[];
for i=1:3
normalized_chen(i,:) = (data_chen(:,i)' - min_chen(i)) / (max_chen(i)-min_chen(i))/2+range_chen(i);
xf1=normalized_chen(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num)==0);
data_Chen=[data_Chen;xf1];
end
%% Chua           
c1=15.6;c2=1;m0=-8/7;m1=-5/7;
c3=33;
  F=@(t,Y)[c1*(Y(2)-Y(1)-(m1*Y(1)+(m0-m1)/2*(abs(Y(1)+1)-abs(Y(1)-1))));...
                 c2*(Y(1)-Y(2)+Y(3));...
                 -c3*Y(2);];  
opt = odeset( 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));     
y0=rand(3,1);
[t,Y] = ode45(F, [0:0.001:100], y0, opt);
data_chua=Y(transient:end,:);
max_chua=zeros(3,1);
min_chua=zeros(3,1);
for i=1:3
max_chua(i) = max(data_chua(:,i));
min_chua(i) = min(data_chua(:,i));
end
normalized_chua=data_chua';
% 归一化
data_Chua=[];
for i=1:3
normalized_chua(i,:) = (data_chua(:,i)' - min_chua(i)) / (max_chua(i)-min_chua(i))/2+range_chua(i);
xf1=normalized_chua(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num)==0);
data_Chua=[data_Chua;xf1];
end
%%

beta_loz=0.1*ones(1,size(data_loz,2));
data_loz=[data_loz;beta_loz];

beta_ros=0.5*ones(1,size(data_loz,2));
data_ros=[data_ros;beta_ros];

% beta_lim=0.2*ones(1,size(data_lim,2));
% data_lim=[data_lim;beta_lim];
% 
% beta_qua=0.2*ones(1,size(data_qua,2));
% data_qua=[data_qua;beta_qua];

beta_Chen=0.9*ones(1,size(data_Chen,2));
data_Chen=[data_Chen;beta_Chen];

beta_Chua=1.3*ones(1,size(data_Chua,2));
data_Chua=[data_Chua;beta_Chua];

%%
train_data=[];
train_data=[data_loz,data_ros,data_Chen,data_Chua];
testdata=[data_loz;data_ros;data_Chen;data_Chua];
% figure
% plot3(train_data(1,:),train_data(2,:),train_data(3,:));