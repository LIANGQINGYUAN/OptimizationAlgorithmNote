%{
���д��ļ���������alpha��ʹ��options(18)=0.5����alpha=0.5
Ȼ������ rs_demo('f_p',[0;-2],options)

��RS�㷨�У���alphaΪ1.5��ʱ���0.5�����״ﵽȫ�����ŵ�
%}

function [x,N]=rs_demo(funcname,xnew,options)
%Naive random search demo
% [x,N]=random_search(funcname,xnew,options);
% print = options(1);
% alpha = options(18);

%nargin=number of input arguments
%�����ж������������
%~=�������ڣ�
if nargin ~= 3
    options = [];
    if nargin ~= 2
        disp('Wrong number of arguments.');
        return;
    end
end
%options(14) ��������ֵ��������(Ĭ��ֵΪ100 ��������)��
if length(options) >= 14
    if options(14)==0
        options(14)=1000*length(xnew);
    end
else
    options(14)=1000*length(xnew);
end
%options(18) ��������(Ĭ��ֵΪ1���С)��
if length(options) < 18
    options(18)=1.0; %optional step size
end

format compact;
format short e;

%foptions ���������Ż��Ŀ��ƣ�matlab���ṩ��18����������Щ�������Ż��Ľ������ߺܹؼ������á�
options = foptions(options);
print = options(1);%options(1) ������ʾ����(Ĭ��ֵΪ0)������1ʱ��ʾһЩ�����
epsilon_x = options(2);%�Ż���x�ľ��ȿ���(Ĭ��ֵΪ1e �C4)��

epsilon_g = options(3);%�Ż�����F�ľ��ȿ���(Ĭ��ֵΪ1e �C4)��
max_iter=options(14);%options(14)��������ֵ��������(Ĭ��ֵΪ100 ��������)��
alpha0 = options(18);%ptions(18)����������(Ĭ��ֵΪ1���С)��
if funcname == 'f_r',
    ros_cnt
elseif funcname == 'f_p',
    pks_cnt;
end %if

%����ʼ��
if length(xnew) == 2
    plot(xnew(1),xnew(2),'o')
    text(xnew(1),xnew(2),'Start Point')
    xlower = [-2;-1];
    xupper = [2;3];
end

%feval���ǰ���֪�����ݻ���Ŵ��뵽һ������õĺ��������
%������Ϊfuncname��Ӧ����������Ϊxnew
f_0=feval(funcname,xnew);
xbestcurr = xnew;
xbestold = xnew;
f_best=feval(funcname,xnew);
f_best=10^(sign(f_best))*f_best;

%��������
for k = 1:max_iter,
    %���㵱ǰ��ĺ���ֵ
    xcurr=xbestcurr;
    f_curr=feval(funcname,xcurr);
    alpha = alpha0;
    % �����µĵ�
    xnew = xcurr + alpha*(2*rand(length(xcurr),1)-1);
    % �ж��µ��Ƿ�Խ��
    for i=1:length(xnew),
        xnew(i) = max(xnew(i),xlower(i));
        xnew(i) = min(xnew(i),xupper(i));
    end %for
    %�жϲ��������ŵ�ͺ���ֵ
    f_new=feval(funcname,xnew);
    if f_new < f_best,
        xbestold = xbestcurr;
        xbestcurr = xnew;
        f_best = f_new;
    end
    
    %��ӡ�������̲���
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
    
    %��ֹ��������
    if norm(xnew-xbestold) <= epsilon_x*norm(xbestold)
        disp('Terminating: Norm of difference between iterates less than');
        disp(epsilon_x);
        break;
    end %if
    
    %��ӡ
    pltpts(xbestcurr,xbestold);
    if k == max_iter
        disp('Terminating with maximum number of iterations');
    end %if

end %for

%nargout�������������������
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