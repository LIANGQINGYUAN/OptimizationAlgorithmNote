options(1)=1;
options(14)=50;
[x,y]=gar([0,10;4,6],20,'f_wave',options);
disp('GA Solution:');
disp(x);
disp('Objective function value:');
disp(y);