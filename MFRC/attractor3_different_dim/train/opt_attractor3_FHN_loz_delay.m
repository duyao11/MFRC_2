%% config
clc
clear
iter_max = 200;
repeat_num =4; % ensemble average size
% 1~5: eig_rho, W_in_a, a, reg, d.
lb = [0 0 0 10^-10 0];
ub = [3 3 1 10^-2  1];
[train_data,testdata,data_num,data_len]= func_dataset_attractor3_delay();
%testdata=func_gen_dynamic_test_data(A_set1,Y0_set);
options = optimoptions('surrogateopt','MaxFunctionEvaluations',iter_max,'PlotFcn','surrogateoptplot');
%options = optimoptions('surrogateopt','MaxFunctionEvaluations',iter_max,'Display','off','PlotFcn',[]);
filename = ['opt_Logi_1_' datestr(now,30) '_' num2str(randi(999)) '.mat'];
min_rmse = @(x) (func_train_repeat_attractor3_delay(x,repeat_num,train_data,testdata,data_num,data_len));
%% main (don't need to change this part)
%rng((now*1000-floor(now*1000))*100000)
tic
[opt_result,opt_fval,opt_exitflag,opt_output,opt_trials] = surrogateopt(min_rmse,lb,ub,options);
toc
save(filename)
if ~ispc
    exit;
end