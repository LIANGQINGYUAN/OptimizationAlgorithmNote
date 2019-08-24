% A particle swarm optimizer
% to find the minimum/maximum of the MATLABs’ peaks function
% D---# of inputs to the function (dimension of problem)
clear
%Parameters
% 参数设置
ps=10;%粒子个数
D=2; %维度
ps_lb=-3;%粒子取值下限
ps_ub=3; %粒子取值上限
vel_lb=-1;%速度取值下限
vel_ub=1; %速度取值上限
iteration_n = 50;%迭代次数
range = [-3, 3; -3, 3]; % 输入值的范围
% 画peaks函数的轮廓
[x, y, z] = peaks;
%pcolor(x,y,z); shading interp; hold on;
%contour(x, y, z, 20, ’r’);
mesh(x,y,z)
% hold off;
%colormap(gray);
set(gca,'Fontsize',14)
axis([-3 3 -3 3 -9 9])
%axis square;
xlabel('x_1','Fontsize',14);
ylabel('x_2','Fontsize',14);
zlabel('f(x_1,x_2)','Fontsize',14);
hold on

upper = zeros(iteration_n, 1);
average = zeros(iteration_n, 1);
lower = zeros(iteration_n, 1);
% initialize population of particles and their velocities at time zero,
% format of pos= (particle#, dimension)
% construct random population positions bounded by VR
% need to bound positions
% 初始化位置和速度
ps_pos=ps_lb + (ps_ub-ps_lb).*rand(ps,D);
% need to bound velocities between -mv,mv
ps_vel=vel_lb + (vel_ub-vel_lb).*rand(ps,D);
% initial pbest positions
p_best = ps_pos;
% returns column of cost values (1 for each particle)
% f1='3*(1-ps_pos(i,1))^2*exp(-ps_pos(i,1)^2-(ps_pos(i,2)+1)^2)';
% f2='-10*(ps_pos(i,1)/5-ps_pos(i,1)^3-ps_pos(i,2)^5)*exp(-ps_pos(i,1)^2-ps_pos(i,2)^2)';
% f3='-(1/3)*exp(-(ps_pos(i,1)+1)^2-ps_pos(i,2)^2)';
% f1+f2+f3为完整的函数

%初始化全局最优和每个粒子经过的最好的位置
p_best_fit=zeros(ps,1);
for i=1:ps
    g1(i)=3*(1-ps_pos(i,1))^2*exp(-ps_pos(i,1)^2-(ps_pos(i,2)+1)^2);
    g2(i)=-10*(ps_pos(i,1)/5-ps_pos(i,1)^3-ps_pos(i,2)^5)*exp(-ps_pos(i,1)^2-ps_pos(i,2)^2);
    g3(i)=-(1/3)*exp(-(ps_pos(i,1)+1)^2-ps_pos(i,2)^2);
    p_best_fit(i)=g1(i)+g2(i)+g3(i);
end
p_best_fit;
hand_p3=plot3(ps_pos(:,1),ps_pos(:,2),p_best_fit','*k','markersize',15);
% initial g_best
[g_best_val,g_best_idx] = max(p_best_fit);
%[g_best_val,g_best_idx] = min(p_best_fit); this is to minimize
g_best=ps_pos(g_best_idx,:);

% get new velocities, positions (this is the heart of the PSO algorithm)
% 迭代步骤
for k=1:iteration_n
    for count=1:ps
        ps_vel(count,:) = 0.729*ps_vel(count,:)... % prev vel
        +1.494*rand*(p_best(count,:)-ps_pos(count,:))... % independent
        +1.494*rand*(g_best-ps_pos(count,:)); % social
    end
    ps_vel;
    % update new position
    ps_pos = ps_pos + ps_vel;
    %update p_best
    for i=1:ps
        g1(i)=3*(1-ps_pos(i,1))^2*exp(-ps_pos(i,1)^2-(ps_pos(i,2)+1)^2);
        g2(i)=-10*(ps_pos(i,1)/5-ps_pos(i,1)^3-ps_pos(i,2)^5)*exp(-ps_pos(i,1)^2-ps_pos(i,2)^2);
        g3(i)=-(1/3)*exp(-(ps_pos(i,1)+1)^2-ps_pos(i,2)^2);
        ps_current_fit(i)=g1(i)+g2(i)+g3(i);
        if ps_current_fit(i)>p_best_fit(i)
            p_best_fit(i)=ps_current_fit(i);
            p_best(i,:)=ps_pos(i,:);
        end
    end
    p_best_fit;
    %update g_best
    [g_best_val,g_best_idx] = max(p_best_fit);
    g_best=ps_pos(g_best_idx,:);
    % Fill objective function vectors
    upper(k) = max(p_best_fit);
    average(k) = mean(p_best_fit);
    lower(k) = min(p_best_fit);
    set(hand_p3,'xdata',ps_pos(:,1),'ydata',ps_pos(:,2),'zdata',ps_current_fit');
    drawnow
    pause%按一下回车键就迭代一次
end

g_best
g_best_val
%画图——最终结果图
figure;
x = 1:iteration_n;
plot(x, upper, 'o', x, average, 'x', x, lower, '*');
hold on;
plot(x, [upper average lower]);
hold off;
legend('Best', 'Average', 'Poorest');
xlabel('Iterations'); ylabel('Objective function value');