clc;clear;
color_f1= [110,69,210]/255;
color_f5=addcolorplus(73);
color_f2= [18,192,190]/255;
color_f3=addcolorplus(13);
color_f7=[218,83,17]/255;
color_f4=addcolorplus(242);
color_f6=addcolorplus(1);


load 4D_ros.mat Y1 testdata
testLen=2000;
%%
figure
plot3(Y1(9,500:testLen),Y1(10,500:testLen),Y1(11,500:testLen),'LineStyle', '-','linewidth',1,'color',color_f3);
hold on;
plot3(testdata(13*1+9,500:testLen),testdata(13*1+10,500:testLen),testdata(13*1+11,500:testLen),'LineStyle', '--','linewidth',1,'color',color_f6);
xlabel('$x$', 'interpreter','latex');
ylabel('$y$','interpreter','latex');
zlabel('$z$', 'interpreter','latex');
%legend('show', 'FontSize', 12, 'Location', 'northwest');
%grid on;
set(gca,'FontName','Times New Roman','FontSize',13);
set(gcf,'unit','centimeters','position',[8 1 8 7]);
left_margin = 0.2;   % 左边距
bottom_margin = 0.15; % 下边距
right_margin = 0.02; % 右边距
top_margin = 0.03;   % 上边距

% 计算坐标轴的位置
width = 1 - left_margin - right_margin;
height = 1 - bottom_margin - top_margin;

set(gca, 'Position', [left_margin, bottom_margin, width, height]);
print('att1.tif', '-dtiff', '-r600');
%%
figure
plot3(Y1(9,500:testLen),Y1(11,500:testLen),Y1(12,500:testLen),'LineStyle', '-','linewidth',1,'color',color_f3);
hold on;
plot3(testdata(13*1+9,500:testLen),testdata(13*1+11,500:testLen),testdata(13*1+12,500:testLen),'LineStyle', '--','linewidth',1,'color',color_f6);
xlabel('$x$', 'interpreter','latex');
ylabel('$z$','interpreter','latex');
zlabel('$w$', 'interpreter','latex');
%legend('show', 'FontSize', 12, 'Location', 'northwest');
%grid on;
set(gca,'FontName','Times New Roman','FontSize',13);
set(gcf,'unit','centimeters','position',[8 1 8 7]);
left_margin = 0.2;   % 左边距
bottom_margin = 0.15; % 下边距
right_margin = 0.02; % 右边距
top_margin = 0.03;   % 上边距

% 计算坐标轴的位置
width = 1 - left_margin - right_margin;
height = 1 - bottom_margin - top_margin;

set(gca, 'Position', [left_margin, bottom_margin, width, height]);
print('att2.tif', '-dtiff', '-r600');
%%
figure
plot3(Y1(10,500:testLen),Y1(11,500:testLen),Y1(12,500:testLen),'LineStyle', '-','linewidth',1,'color',color_f3);
hold on;
plot3(testdata(13*1+10,500:testLen),testdata(13*1+11,500:testLen),testdata(13*1+12,500:testLen),'LineStyle', '--','linewidth',1,'color',color_f6);
xlabel('$y$', 'interpreter','latex');
ylabel('$z$','interpreter','latex');
zlabel('$w$', 'interpreter','latex');
%legend('show', 'FontSize', 12, 'Location', 'northwest');
%grid on;
set(gca,'FontName','Times New Roman','FontSize',13);
set(gcf,'unit','centimeters','position',[8 1 8 7]);
left_margin = 0.2;   % 左边距
bottom_margin = 0.15; % 下边距
right_margin = 0.02; % 右边距
top_margin = 0.03;   % 上边距

% 计算坐标轴的位置
width = 1 - left_margin - right_margin;
height = 1 - bottom_margin - top_margin;

set(gca, 'Position', [left_margin, bottom_margin, width, height]);
print('att3.tif', '-dtiff', '-r600');
%%
load FHN.mat Y1 testdata
figure
%hold on;
plot(Y1(11,1000:testLen),Y1(12,1000:testLen),'LineStyle', '-','linewidth',1,'color',color_f3);
 hold on;
plot(testdata(13*0+11,200:testLen),testdata(13*0+12,200:testLen),'LineStyle', '--','linewidth',1,'color',color_f6);
 set(gca,'FontName','Times New Roman','FontSize',13);
 xlabel('$x$', 'interpreter','latex');
ylabel('$y$','interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13);
set(gcf,'unit','centimeters','position',[8 1 8 7]);
left_margin = 0.18;   % 左边距
bottom_margin = 0.2; % 下边距
right_margin = 0.02; % 右边距
top_margin = 0.03;   % 上边距

% 计算坐标轴的位置
width = 1 - left_margin - right_margin;
height = 1 - bottom_margin - top_margin;

set(gca, 'Position', [left_margin, bottom_margin, width, height]);
xlim([0.19, 0.71]); % Adjust to match the Python plot more closely
ylim([0.19, 0.73]);
print('att4.tif', '-dtiff', '-r600');
 %%%
load loz.mat Y1 testdata
testLen=2000;
figure
%hold on;
plot(Y1(11,200:testLen),Y1(12,200:testLen),'LineStyle', '-','linewidth',1,'color',color_f3);
 hold on;
plot(testdata(13*2+11,200:testLen),testdata(13*2+12,200:testLen),'LineStyle', '--','linewidth',1,'color',color_f6);
 set(gca,'FontName','Times New Roman','FontSize',13);
 xlabel('$y$', 'interpreter','latex');
ylabel('$z$','interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13);
set(gcf,'unit','centimeters','position',[8 1 8 7]);
left_margin = 0.18;   % 左边距
bottom_margin = 0.2; % 下边距
right_margin = 0.02; % 右边距
top_margin = 0.03;   % 上边距

% 计算坐标轴的位置
width = 1 - left_margin - right_margin;
height = 1 - bottom_margin - top_margin;

set(gca, 'Position', [left_margin, bottom_margin, width, height]);
xlim([-0.2, 0.3]); % Adjust to match the Python plot more closely
ylim([-0.2, 0.3]);
print('att5.tif', '-dtiff', '-r600');
%%
figure
%hold on;
plot(Y1(10,200:testLen),Y1(11,200:testLen),'LineStyle', '-','linewidth',1,'color',color_f3);
hold on;
plot(testdata(13*2+10,200:testLen),testdata(13*2+11,200:testLen),'LineStyle', '--','linewidth',1,'color',color_f6);
 set(gca,'FontName','Times New Roman','FontSize',13);
 xlabel('$x$', 'interpreter','latex');
ylabel('$y$','interpreter','latex');
set(gca,'FontName','Times New Roman','FontSize',13);
set(gcf,'unit','centimeters','position',[8 1 8 7]);
left_margin = 0.18;   % 左边距
bottom_margin = 0.2; % 下边距
right_margin = 0.02; % 右边距
top_margin = 0.03;   % 上边距

% 计算坐标轴的位置
width = 1 - left_margin - right_margin;
height = 1 - bottom_margin - top_margin;

set(gca, 'Position', [left_margin, bottom_margin, width, height]);
xlim([-0.2, 0.3]); % Adjust to match the Python plot more closely
ylim([-0.2, 0.3]);
print('att6.tif', '-dtiff', '-r600');