function dec = bin2dec(bin,range);
%function dec = bin2dec(bin,range);
%Function to convert from binary (bin) to decimal (dec) in a given range
index = polyval(bin,2);
dec = index*((range(2)-range(1))/(2^length(bin)-1)) + range(1);