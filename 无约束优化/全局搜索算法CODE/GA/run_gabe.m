clear;
options(1)=1;
[x,y]=ga_be(8,10,'fit_func1',options);
f='f_manymax';
range=[-10,10];
disp('GA Solution:');
disp(bin2dec(x,range));
disp('Objective function value:');
disp(y);

