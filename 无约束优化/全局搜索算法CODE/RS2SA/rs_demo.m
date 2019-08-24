%{
运行此文件首先设置alpha，使用options(18)=0.5设置alpha=0.5
然后运行 rs_demo('f_p',[0;-2],options)

在RS算法中，当alpha为1.5的时候比0.5更容易达到全局最优点
%}

function [x,N]=rs_demo(funcname,xnew,options)
%Naive random search demo
% [x,N]=random_search(funcname,xnew,options);
% print = options(1);
% alpha = options(18);

%nargin=number of input arguments
%用于判断输入参数个数
%~=（不等于）
if nargin ~= 3
    options = [];
    if nargin ~= 2
        disp('Wrong number of arguments.');
        return;
    end
end
%options(14) ：函数求值的最大次数(默认值为100 变量个数)。
if length(options) >= 14
    if options(14)==0
        options(14)=1000*length(xnew);
    end
else
    options(14)=1000*length(xnew);
end
%options(18) 步长设置(默认值为1或更小)。
if length(options) < 18
    options(18)=1.0; %optional step size
end

format compact;
format short e;

%foptions 函数对于优化的控制，matlab共提供了18个参数，这些参数对优化的进行起者很关键的作用。
options = foptions(options);
print = options(1);%options(1) 参数显示控制(默认值为0)。等于1时显示一些结果。
epsilon_x = options(2);%优化点x的精度控制(默认值为1e –4)。

epsilon_g = options(3);%优化函数F的精度控制(默认值为1e –4)。
max_iter=options(14);%options(14)：函数求值的最大次数(默认值为100 变量个数)。
alpha0 = options(18);%ptions(18)：步长设置(默认值为1或更小)。
if funcname == 'f_r',
    ros_cnt
elseif funcname == 'f_p',
    pks_cnt;
end %if

%画初始点
if length(xnew) == 2
    plot(xnew(1),xnew(2),'o')
    text(xnew(1),xnew(2),'Start Point')
    xlower = [-2;-1];
    xupper = [2;3];
end

%feval就是把已知的数据或符号带入到一个定义好的函数句柄中
%即给名为funcname相应参数，参数为xnew
f_0=feval(funcname,xnew);
xbestcurr = xnew;
xbestold = xnew;
f_best=feval(funcname,xnew);
f_best=10^(sign(f_best))*f_best;

%迭代过程
for k = 1:max_iter,
    %计算当前点的函数值
    xcurr=xbestcurr;
    f_curr=feval(funcname,xcurr);
    alpha = alpha0;
    % 产生新的点
    xnew = xcurr + alpha*(2*rand(length(xcurr),1)-1);
    % 判断新点是否越界
    for i=1:length(xnew),
        xnew(i) = max(xnew(i),xlower(i));
        xnew(i) = min(xnew(i),xupper(i));
    end %for
    %判断并更新最优点和函数值
    f_new=feval(funcname,xnew);
    if f_new < f_best,
        xbestold = xbestcurr;
        xbestcurr = xnew;
        f_best = f_new;
    end
    
    %打印迭代过程参数
    if print,
        disp('Iteration number k =')
        disp(k); %print iteration index k
        disp('alpha =');
        disp(alpha); %print alpha
        disp('New point =');
        disp(xnew'); %print new point
        disp('Function value =');
        disp(f_new); %print func value at new point
    end %if
    
    %终止条件设置
    if norm(xnew-xbestold) <= epsilon_x*norm(xbestold)
        disp('Terminating: Norm of difference between iterates less than');
        disp(epsilon_x);
        break;
    end %if
    
    %打印
    pltpts(xbestcurr,xbestold);
    if k == max_iter
        disp('Terminating with maximum number of iterations');
    end %if

end %for

%nargout：函数输出参数的数量
if nargout >= 1
    x=xnew;
    if nargout == 2
        N=k;
    end
else
    disp('Final point =');
    disp(xbestcurr');
    disp('Number of iterations =');
    disp(k);
end %if