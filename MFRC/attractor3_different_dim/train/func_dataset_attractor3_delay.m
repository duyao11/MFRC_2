%clc;clear;
function [train_data,testdata,data_num,data_len]= func_dataset_attractor3_delay()
data_num=3;
data_len=ones(1,data_num)*3000;
range=0.2;
%%
range_FHN=[0.2,0.2];
range_4D_ros=[0.2,-0.2,0.2,-0.2];
range_loz=[-0.2,-0.2,-0.2];
%%
all_dim=12;
dim_FHN=all_dim/2;
dim_4D_ros=all_dim/4;
dim_loz=all_dim/3;
mod_num_FHN=200;
mod_num_4D_ros=100;
mod_num_loz=20;
%%
transient=80000;
%% FHN
    epsilon = 1; % 小的时间尺度参数
    a = 0.7;        % 参数a，影响nullcline
    b = 0.4;        % 参数b，影响nullcline
    I = 0.5;        % 外部输入电流
    % 初始条件
    v0 = -1;        % 初始膜电位
    w0 = 1;         % 初始恢复变量
    y0 = [v0; w0];
    % 时间范围
    tspan = [0:0.001: 1000];

    % 使用ode45求解ODE
    [t, Y] = ode45(@(t, y) fhnODE(t, y, epsilon, a, b, I), tspan, y0);

data_FHN=Y(transient:end,:);
max_FHN=zeros(2,1);
min_FHN=zeros(2,1);
for i=1:2
max_FHN(i) = max(data_FHN(:,i));
min_FHN(i) = min(data_FHN(:,i));
end
normalized_FHN=data_FHN';
% 归一化
Data_FHN=[];
for i=1:2
normalized_FHN(i,:) = (data_FHN(:,i)' - min_FHN(i)) / (max_FHN(i)-min_FHN(i))/2+range;
xf1=normalized_FHN(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num_FHN)==0);
Data_FHN=[Data_FHN;xf1];
end
Dataset_FHN=zeros(6,length(xf1)-(dim_FHN-1));
for i=1:dim_FHN
    Dataset_FHN(2*i-1,:)=Data_FHN(1,i:length(xf1)-(dim_FHN-i));
    Dataset_FHN(2*i,:)=Data_FHN(2,i:length(xf1)-(dim_FHN-i));
end
%% 4D ros
transient=80003;
tspan = [0:0.001:1000];
% Initial conditions
x0 = [0, 1, 0, 1];  % Initial conditions for [x, y, z, w]
[t, X] = ode45(@hyperchaotic_system, tspan, x0);
data_4D_rosser=X(transient:end,:);
max_4D_ros=zeros(4,1);
min_4D_ros=zeros(4,1);
for i=1:4
max_4D_ros(i) = max(data_4D_rosser(:,i));
min_4D_ros(i) = min(data_4D_rosser(:,i));
end
normalized_4D_ros=data_4D_rosser';
% 归一化
Data_4D_ros=[];
for i=1:4
normalized_4D_ros(i,:) = (data_4D_rosser(:,i)' - min_4D_ros(i)) / (max_4D_ros(i)-min_4D_ros(i))/2-range_4D_ros(i);
xf1=normalized_4D_ros(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num_4D_ros)==0);
Data_4D_ros=[Data_4D_ros;xf1];
end
Dataset_4D_ros=zeros(all_dim,length(xf1)-(dim_4D_ros-1));
for i=1:dim_4D_ros
    Dataset_4D_ros(4*i-3,:)=Data_4D_ros(1,i:length(xf1)-(dim_4D_ros-i));
    Dataset_4D_ros(4*i-2,:)=Data_4D_ros(2,i:length(xf1)-(dim_4D_ros-i));
    Dataset_4D_ros(4*i-1,:)=Data_4D_ros(3,i:length(xf1)-(dim_4D_ros-i));
    Dataset_4D_ros(4*i,:)=Data_4D_ros(4,i:length(xf1)-(dim_4D_ros-i));
end
%% Lorenz
transient=80003;
F=@(t,Y)[10*Y(2)-10*Y(1);...
         -Y(1).*Y(3)+28*Y(1)-Y(2);...
         Y(1).*Y(2)-8*Y(3)/3;...
         ];
         
opt = odeset( 'RelTol', 10.0^(-7), 'AbsTol' , 10.0^(-7));     
y0=rand(3,1);
[t,Y] = ode45(F, [0:0.001:500], y0, opt);
data_lorenz=Y(transient:end,:);
max_loz=zeros(3,1);
min_loz=zeros(3,1);
for i=1:3
max_loz(i) = max(data_lorenz(:,i));
min_loz(i) = min(data_lorenz(:,i));
end
normalized_loz=data_lorenz';
% 归一化
Data_loz=[];
for i=1:3
normalized_loz(i,:) = (data_lorenz(:,i)' - min_loz(i)) / (max_loz(i)-min_loz(i))/2-range;
xf1=normalized_loz(i,:);
x_pos=1:size(xf1,2);
xf1=xf1(mod(x_pos,mod_num_loz)==0);
Data_loz=[Data_loz;xf1];
end
Dataset_loz=zeros(all_dim,length(xf1)-(dim_loz-1));
for i=1:dim_loz
    Dataset_loz(3*i-2,:)=Data_loz(1,i:length(xf1)-(dim_loz-i));
    Dataset_loz(3*i-1,:)=Data_loz(2,i:length(xf1)-(dim_loz-i));
    Dataset_loz(3*i,:)=Data_loz(3,i:length(xf1)-(dim_loz-i));
end
%%
beta_FHN=0.5*ones(1,size(Dataset_FHN,2));
Dataset_FHN=[Dataset_FHN;beta_FHN];

beta_4D_ros=1*ones(1,size(Dataset_4D_ros,2));
Dataset_4D_ros=[Dataset_4D_ros;beta_4D_ros];

beta_loz=1.5*ones(1,size(Dataset_loz,2));
Dataset_loz=[Dataset_loz;beta_loz];

%%
train_data=[];
train_data=[Dataset_FHN(:,1:data_len(1)),Dataset_4D_ros(:,1:data_len(1)),Dataset_loz(:,1:data_len(1))];
testdata=[Dataset_FHN(:,1:data_len(1));Dataset_4D_ros(:,1:data_len(1));Dataset_loz(:,1:data_len(1))];
figure
plot(train_data(3,:));
figure
plot3(train_data(1,:),train_data(2,:),train_data(3,:));
function dy = fhnODE(t, y, epsilon, a, b, I)
    v = y(1);
    w = y(2);
    dvdt = v - (v^3)/3 - w + I;       % 膜电位的变化率
    dwdt = epsilon * (v + a - b * w); % 恢复变量的变化率
    dy = [dvdt; dwdt];
end
function dxdt = hyperchaotic_system(t, x)
    a = 0.1;
    b = 0.1;
    c = 14;
    pz = -20;
    pw = -0.28;
    
    dxdt = zeros(4,1);
    dxdt(1) = -x(2) - x(3) + x(4);                           % dx/dt
    dxdt(2) = x(1) + a*x(2) - x(3);                          % dy/dt
    dxdt(3) = b + x(3)*(x(1) - c) + pz*x(1)*x(4);            % dz/dt
    dxdt(4) = -x(1)*x(2) + pw*x(4);                          % dw/dt
end
end