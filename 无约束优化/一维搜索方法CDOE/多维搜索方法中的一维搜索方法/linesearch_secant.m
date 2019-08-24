function alpha=linesearch_secant(grad,x,d)
%Line search using secant method
%{
gradΪ����
xΪ��ʼ��
dΪ������������

ʹ��Ŀ�꺯�����ＫС��Ϊһά������Ŀ��
alphaΪһά���������ͷ���ֵ�������ڵ���x���������alpha��
%}
epsilon=10^(-4); %line search tolerance
max = 100; %maximum number of iterations
alpha_curr=0;%��ǰalphaֵ
alpha=0.001; %��һ��alphaֵ
dphi_zero=feval(grad,x)'*d; % ��ʼ����ֵ�ĵ�������������ĳ˻�
dphi_curr=dphi_zero;        

i=0;
while abs(dphi_curr)>epsilon*abs(dphi_zero),
    %alphaֵ����
    alpha_old=alpha_curr;
    alpha_curr=alpha;
    %����ֵ����
    dphi_old=dphi_curr;
    dphi_curr=feval(grad,x+alpha_curr*d)'*d;
    %���߷�
    alpha=(dphi_curr*alpha_old-dphi_old*alpha_curr)/(dphi_curr-dphi_old);
    %��ֹ�����ж�
    i=i+1;
    if (i >= max) & (abs(dphi_curr)>epsilon*abs(dphi_zero)),
        disp('Line search terminating with number of iterations:');
        disp(i);
        break;
    end
end %while