%{
在两个初始点作为输入的情况下这个函数具有两个局部最小点，这个情况取决于lambda的取值。
而lambda的取值也决定了初始的单纯形，通过使用不同的lambda值，可以从相同的初始点开始达到两个最小值。
在下面的解决方案中，初始的单纯形比较小，因为lambda 0.1。
%}
function [ output_args ] = nm_simplex( input_args )
%Nelder-Mead simplex method
%Based on the program by the Spring 2007 ECE580 student, Hengzhou Ding
disp ('We minimize a function using the Nelder-Mead method.')
disp ('There are two initial conditions.')
disp ('You can enter your own starting point.')
disp ('---------------------------------------------')

% disp(’Select one of the starting points’)
% disp (’[0.55;0.7] or [-0.9;-0.5]’)
% x0=input(’’)
disp (' ')
clear
close all;

disp('Select one of the starting points, or enter your own point')
disp('[0.55;0.7] or [-0.9;-0.5]')
disp('(Copy one of the above points and paste it at the prompt)')
x0=input(' ');

hold on
axis square
%Plot the contours of the objective function
%画目标函数的轮廓
[X1,X2]=meshgrid(-1:0.01:1);
Y=(X2-X1).^4+12.*X1.*X2-X1+X2-3;
[C,h] = contour(X1,X2,Y,20);
clabel(C,h);


% Initialize all parameters
%初始化参数
%{
lambda用于生成n+1个点
rho用于生成pr
chi用于生成pe
gamma用于生成pc
sigma用于生成重点
e1、e2为空间的标准基
%}
lambda=0.1;
rho=1;
chi=2;
gamma=1/2;
sigma=1/2;
e1=[1 0]';
e2=[0 1]';
%x0=[0.55 0.7]’;
%x0=[-0.9 -0.5]’;

% Plot initial point and initialize the simplex
%画初始的迭代点
%x(:,3)=x0表示矩阵x的第三列为向量x0
plot(x0(1),x0(2),'--*');
x(:,3)=x0;
x(:,1)=x0+lambda*e1;
x(:,2)=x0+lambda*e2;

%{
递归迭代的过程
%}
while 1
    % Check the size of simplex for stopping criterion
    % 检查单纯形是否满足停止条件
    simpsize=norm(x(:,1)-x(:,2))+norm(x(:,2)-x(:,3))+norm(x(:,3)-x(:,1));
    if(simpsize<1e-6)
        break;
    end
    %上一次的迭代点
    lastpt=x(:,3);
    % Sort the simplex
    % 对单纯形中的点进行排序——排序后x(:,3)中的点对应函数值最大
    x=sort_points(x,3);
    % Reflection
    % 反射操作
    centro=1/2*(x(:,1)+x(:,2));
    xr=centro+rho*(centro-x(:,3));
    % Accept condition
    % 接受条件判断
    % 条件1：反射点pr对应函数值位于fnl和fs之间
    if(obj_fun(xr)>=obj_fun(x(:,1)) && obj_fun(xr)<obj_fun(x(:,2)))
        x(:,3)=xr;
        % Expand condition
    %条件2：反射点pr对应的函数值小于fs
    elseif(obj_fun(xr)<obj_fun(x(:,1)))
        xe=centro+rho*chi*(centro-x(:,3));
            if(obj_fun(xe)<obj_fun(xr))
                x(:,3)=xe;
            else
                x(:,3)=xr;
            end

    % Outside contraction or shrink
    %条件3：反射点pr对应的函数值大于 fnl，但小于fl。——外收缩
    elseif(obj_fun(xr)>=obj_fun(x(:,2)) && obj_fun(xr)<obj_fun(x(:,3)))
        xc=centro+gamma*rho*(centro-x(:,3));
            if(obj_fun(xc)<obj_fun(x(:,3)))
                x(:,3)=xc;
            else
                x=shrink(x,sigma);
            end
        % Inside contraction or shrink
    % 条件4：反射点pr对应的函数值大于 fl——内收缩
    else
        xcc=centro-gamma*(centro-x(:,3));
        if(obj_fun(xcc)<obj_fun(x(:,3)))
            x(:,3)=xcc;
        else
            x=shrink(x,sigma);
        end
    end
    % Plot the new point and connect
    % 打印反射出来的新的迭代点，并与上一个点连线
    plot([lastpt(1),x(1,3)],[lastpt(2),x(2,3)],'--*');
end
% Output the final simplex (minimizer)
% 输出最终的单纯形
x(:,1)


% obj_fun
%目标函数
function y = obj_fun(x)
y=(x(1)-x(2))^4+12*x(1)*x(2)-x(1)+x(2)-3;

% sort_points
% 单纯形中的点进行排序——最终x(:,1)>x(:,2)>x(:,3)
function y = sort_points(x,N)
for i=1:(N-1)
    for j=1:(N-i)
        if(obj_fun(x(:,j))>obj_fun(x(:,j+1)))
            tmp=x(:,j);
            x(:,j)=x(:,j+1);
            x(:,j+1)=tmp;
        end
    end
end
y=x;

% shrink
% 收缩函数
function y = shrink(x,sigma)
x(:,2)=x(:,1)+sigma*(x(:,2)-x(:,1));
x(:,3)=x(:,1)+sigma*(x(:,3)-x(:,1));
y=x;
