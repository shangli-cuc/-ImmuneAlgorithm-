%一个旅行商人要拜访31个城市，需要选择最短的路径
%种群数目NP=200
%免疫个体维数N=31
%迭代次数G=1000
%克隆个体个数为Ncl=10
%任意两个城市的距离矩阵D

%免疫算法解决TSP问题
%初始化
clear all; %清除所有变量
close all; %清图
clc ;      %清屏
%设置了31个城市的坐标，31*2的矩阵
C=[
    1304 2312;
    3639 1315;
    4177 2244;
    3712 1399;
    3488 1535;
    3326 1556;
    3238 1229;
    4196 1004;
    4312 790;
    4386 570;
    3007 1970;
    2562 1756;
    2788 1491;
    2381 1676;
    1332 695;
    3715 1678;
    3918 2179;
    4061 2370;
    3780 2212;
    3676 2578;
    4029 2838;
    4263 2931;
    3429 1908;
    3507 2367;
    3394 2643;
    3439 3201;
    2935 3240;
    3140 3550;
    2545 2357;
    2778 2826;
    2370 2975
    ];
N=size(C,1); %TSP问题的规模，即城市数目N=31
D=zeros(N); %任意两个城市距离间隔矩阵,初始化都为0  N*N的零矩阵

%任意两个城市距离间隔矩阵
for i=1:N
    for j=1:N
        D(i,j)=((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;  %第i行第j列表示第i个城市与第j个城市的距离
    end
end

NP=200; %免疫个体数目 （初代解的个数）
G=1000; %最大免疫代数（迭代次数）
f=zeros(N,NP); %用于存储初始抗体，每一列都是一个解向量  31*200矩阵

for i=1:NP
    f(:,i)=randperm(N); %随机生成初始抗体  f矩阵的第i列 31*200矩阵
end

len=zeros(NP,1);  %存储路径长度，解的能力，路径长度越短，能力越强  200*1零矩阵
for i=1:NP
    len(i)=func_length(D,f(:,i),N); %随机生成初始种群  第i种解的路径长度
end

[Sortlen,Index]=sort(len);  %Sortlen是对200种解的结果由小到大排序Index是索引数组，size(Index)==size(len)
SortFirst=f(:,Index); %种群个体排序 31*200矩阵内容和f矩阵一样，但是总路径len越短的解越靠前排列
generation=0;               %免疫代数(迭代次数)
Nc1=10;             %克隆个数

%%%%%%%%%%%%%%%%%%%%%%免疫循环开始%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while generation<G
    for i=1:NP/2    %选激励度前NP/2个个体进行免疫操作，激励度指抗体应答抗原的综合能力，与亲和度成正比，这里是选择总路径len较短的前100种解
        
        a=SortFirst(:,i);%f的前NP/2列随机产生的第i个解 31*1矩阵
        Ca=repmat(a,1,Nc1);%将第i个解克隆Nc1，即10个，组成31*10矩阵
        for j=1:Nc1
            p1=floor(1+N*rand());%rand随机产生1个0到1之间的数，floor向下取整，p1，p2在1~31之间
            p2=floor(1+N*rand());
            while p1==p2
                p1=floor(1+N*rand());
                p2=floor(1+N*rand());
            end
            %当两个随机数不等时，进行元素的交换，改变一个解中某两个城市的次序，相当于原来相同的10列解进行了10次不同的交叉变异 Cross and Aberrance
            temp=Ca(p1,j);
            Ca(p1,j)=Ca(p2,j);
            Ca(p2,j)=temp;
        end
        Ca(:,1)=SortFirst(:,i);   %保留克隆源个体 将交叉变异前Sortf的解赋值给Ca的第1列，现在每种解只剩下9次不同的交叉变异，总共有9*10种交叉变异
        %%%%%%%%%%%%%%%%克隆抑制，保留亲和度最高的个体%%%%%%%%%%%%%%%%
        for j=1:Nc1
            Calen(j)=func_length(D,Ca(:,j),N);%交叉变异后解的路径长度Calen，第1列是原本未变异的解的路径长度，其余9列是变异后的路径长度
        end
        [SortCalen,Index]=sort(Calen);%由小到大排序
        SortCa=Ca(:,Index);%将小的放在前面，大的放在后面
        af(:,i)=SortCa(:,1);%将路径最短的解挑选出来
        alen(i)=SortCalen(1);%记录这种解的最短路径
    end%最后保留100种解存储在af中 31*100矩阵
    %%%%%%%%%%%%%%%%%种群刷新%%%%%%%%%%%
    for i=1:NP/2
        bf(:,i)=randperm(N); %随机生成初始种群 31*100矩阵
        blen(i)=func_length(D,bf(:,i),N); %计算路径长度 100*1矩阵
    end
    %%%%%%%%%%%%%%%%%免疫种群与新种群合并%%%%%%%%%%%
    f=[af,bf];%构成新的解的集合 31*200
    len=[alen,blen];
    [Sortlen,Index]=sort(len);%再次排序
    SortFirst=f(:,Index);%排序后的解集合
    generation=generation+1;%迭代次数加一
    trace(generation)=Sortlen(1);%将最小路径赋值给trace
end %每次迭代交叉变异900次，迭代1000次

%%%%%%%%%%%%%%%%%%%%%%输出优化结果%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BestResult=SortFirst(:,1);   %最优变量 31*1
BestLen=trace(end);  %最优值
%图1
figure
for i=1:N-1
    plot([C(BestResult(i),1),C(BestResult(i+1),1)],[C(BestResult(i),2),C(BestResult(i+1),2)],'bo-');
    hold on;
end

plot([C(BestResult(N),1),C(BestResult(1),1)],[C(BestResult(N),2),C(BestResult(1),2)],'ro-');
title(['优化最短距离:',num2str(BestLen)]);
%图2亲和度曲线，用每次迭代得到的最优解作为亲和度
figure
% for i=1:1000
%     trace(i)=log(trace(i));
% % end
% x=1:1000;
% plot(trace)
for i=1:generation-1
    plot([i,i+1],[log(trace(i)),log(trace(i+1))],'bo-');
    hold on;
end
xlabel('迭代次数')
ylabel('目标函数值')
title('亲和度进化曲线')




