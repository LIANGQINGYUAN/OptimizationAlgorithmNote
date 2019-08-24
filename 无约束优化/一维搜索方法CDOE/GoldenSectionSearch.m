%Matlab routine for Golden Section Search
left=1;
right=2;
uncert=0.23;%最终的压缩长度
rho=(3-sqrt(5))/2;

N=ceil(log(uncert/(right-left))/log(1-rho)) %print N

lower='a';
a=left+(1-rho)*(right-left);
f_a=f(a);

for i=1:N,
    if lower=='a'
        b=a
        f_b=f_a
        a=left+rho*(right-left)
        f_a=f(a)
    else
        a=b
        f_a=f_b
        b=left+(1-rho)*(right-left)
        f_b=f(b)
    end %if
    if f_a<f_b
        right=b;
        lower='a'
    else
        left=a;
        lower='b'
    end %if
    New_Interval = [left,right]
end %for i