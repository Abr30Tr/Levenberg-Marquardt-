%%%%%%%%%%%% how to use %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to use the class user should provide initial value (x0),provide it as a row vector 
% an anonymous function (fun) of residual equation(least square problem) and its Jacobian (jac). 
% prepare the data to be fitted in single array to provide as an argument
% pass the above element as an argument to solve_inc and save_result functions 
% save_result provide the update of X, use the end iterative value of array to get
% the best fit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('DataFilteredlong.txt');
f = @(t, x1, x2) x1*t + x2;
fun = @ (y,t,x) y-(x(1)*t + x(2));
jac = @(x,t) -[t, ones(length(t),1)];
%--------------------------------------------------------------------------
% fitting part one of the data 
%--------------------------------------------------------------------------
% the first time interval values of time and voltage for linear fitting 
time1 = DataFilteredlong(1:242,1);
volt1 = DataFilteredlong(1:242,2);
% %--------------------------------------------------------------------------
t= time1;
y = volt1;
x0       = [0.1622 0.7943];
%--------------------------------------------------------------------------
fitData1(:,1) = t;
fitData1(:,2) = y;
x = x0(:);  
LQ = levemberg_marq();
LQ.set_tol_esp1( 1e-3 );
[sol r_sol A_res G_res meu norm_G] = LQ.solve_inc(fun,jac,x,fitData1);
[x lm_result  hist] = LQ.save_result(fun,jac,sol,r_sol,A_res,G_res,norm_G,meu,x,fitData1);

plot(t,y,'r','LineWidth',1);hold on
hold all
%--------------------------------------------------------------------------
% plot for the linear part
%--------------------------------------------------------------------------
y_sol = f(t, x(1,end),x(2,end));
s_1 = f(t(end), x(1,end),x(2,end));
plot(t, y_sol,'b','LineWidth',1); hold on
 grid;
legend('sample data', 'fitted function', 'Location','East')
%--------------------------------------------------------------------------
% fitting part two of the data 
%--------------------------------------------------------------------------
clear fitData
%--------------------------------------------------------------------------
% the second time interval values of time and voltage for linear fitting 
time2 = DataFilteredlong(242:528,1);
volt2 = DataFilteredlong(242:528,2);
%--------------------------------------------------------------------------
t_1= time2;
y_1 = volt2;
x_0       = [-100];
fitData2(:,1) = t_1;
fitData2(:,2) = y_1;
x_1 = x_0(:);
f_1 = @(t_1,x1) 1.189652719636752e+02 - x1*0.9942*10e-4 + x1*t_1;
fun = @(y_1,t_1,x_1) y_1 - 1.189652719636752e+02 + x_1(1)*0.9942*10e-4 - x_1(1)*t_1;
jac = @(x_1,t_1) [4971/5000000 - t_1];
%--------------------------------------------------------------------------
[sol r_sol A_res G_res meu norm_G] = LQ.solve_inc(fun,jac,x_1,fitData2);
% fprintf('x = %g\n',x);
[x_1 lm_result  hist] = LQ.save_result(fun,jac,sol,r_sol,A_res,G_res,norm_G,meu,x_1,fitData2);
%--------------------------------------------------------------------------
plot(t_1,y_1,'r','LineWidth',1);hold on
hold all
%--------------------------------------------------------------------------
% fitting part three of the data 
%--------------------------------------------------------------------------
y_1 = f_1(t_1, x_1(1,end));
plot(t_1,y_1,'b','LineWidth',1); hold on
s_2 = -(t_1(1)*x_1(1,end)) + t_1(end)*x_1(1,end);
s_3 = log(s_1 + s_2);
grid;
legend('sample data', 'fitted function', 'Location','East')
clear fitData
%--------------------------------------------------------------------------
f_3 =  @(t_3,x1) exp(x1.*t_3-x1.*t_3(1)+4.040748722489156);
fun =  @(y,t_3,x_3) y-exp(x_3(1).*t_3-x_3(1).*t_3(1)+4.040748722489156);
jac =  @(x_3,t_3)[-exp(x_3.*t_3 - (1212181669943639.*x_3)/1152921504606846976 + 4549478610224993/1125899906842624).*(t_3 - 1212181669943639/1152921504606846976)];
%--------------------------------------------------------------------------
% the third time interval values of time and voltage for nonlinear fitting 
time7 = DataFilteredlong(528:end,1);
volt7 = DataFilteredlong(528:end,2);
%--------------------------------------------------------------------------
t_3= time7;
y_3 = volt7;
x_0_3 = [1800000];
fitData3(:,1) = t_3;
fitData3(:,2) = y_3;
x_3 = x_0_3(:); 
%--------------------------------------------------------------------------
[sol r_sol A_res G_res meu norm_G] = LQ.solve_inc(fun,jac,x_3,fitData3);
[x_3 lm_result  hist] = LQ.save_result(fun,jac,sol,r_sol,A_res,G_res,norm_G,meu,x_3,fitData3);
%--------------------------------------------------------------------------
plot(t_3,y_3,'r','LineWidth',1);hold on
hold all
%--------------------------------------------------------------------------
y_3 = f_3(t_3,x_3(1,end));
plot(t_3,y_3,'b','LineWidth',1); hold on
grid;
legend('sample data', 'fitted function', 'Location','East')
%--------------------------------------------------------------------------