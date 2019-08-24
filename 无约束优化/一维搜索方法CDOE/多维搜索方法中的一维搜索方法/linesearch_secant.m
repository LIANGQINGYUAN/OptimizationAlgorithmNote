function alpha=linesearch_secant(grad,x,d)
%Line search using secant method
%{
grad为函数
x为初始点
d为方向搜索方向

使得目标函数到达极小点为一维搜索的目标
alpha为一维搜索参数和返回值（类似于迭代x，这里迭代alpha）
%}
epsilon=10^(-4); %line search tolerance
max = 100; %maximum number of iterations
alpha_curr=0;%当前alpha值
alpha=0.001; %下一个alpha值
dphi_zero=feval(grad,x)'*d; % 初始函数值的导数与搜索方向的乘积
dphi_curr=dphi_zero;        

i=0;
while abs(dphi_curr)>epsilon*abs(dphi_zero),
    %alpha值设置
    alpha_old=alpha_curr;
    alpha_curr=alpha;
    %函数值设置
    dphi_old=dphi_curr;
    dphi_curr=feval(grad,x+alpha_curr*d)'*d;
    %割线法
    alpha=(dphi_curr*alpha_old-dphi_old*alpha_curr)/(dphi_curr-dphi_old);
    %终止条件判断
    i=i+1;
    if (i >= max) & (abs(dphi_curr)>epsilon*abs(dphi_zero)),
        disp('Line search terminating with number of iterations:');
        disp(i);
        break;
    end
end %while