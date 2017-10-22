function [ len ] = func_length( D,f,N )
len=D(f(N),f(1));%f(N)表示f矩阵的从上往下，从左往右数第N个元素，这里是第31个城市到第1个城市的距离赋给len
for i=1:N-1
    len=len+D(f(i),f(i+1)); %按照随机分配的31个城市顺序计算总的路径长度
end
end

