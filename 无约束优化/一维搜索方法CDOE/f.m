% f.m
function y=f(x)
y=8*exp(1-x)+7*log(x);

%{
����-��ͼ
fplot('f',[1 2]);
xlabel('x');
ylabel('f(x)');
%}