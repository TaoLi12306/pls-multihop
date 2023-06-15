function u= F_two_hop_independent(Rs,Rv,Am,i,lambdab,lambdae)
%F_TWO_HOP 此处显示有关此函数的摘要
%   此处显示详细说明
pm=zipf(10,i);%内容m的请求概率
pt=1000;     %发射功率
sigma=1 ; %噪声功率
sigma2=1;
a=4;       %衰减指数
Rs=Rs;      %传输速率
lambdab=lambdab;  %基站密度
lambdae=lambdae;    %窃听者的密度
Rv=Rv;        %冗余速率
am=Am; %存储概率 


sum1=0;
sum2=0;

d1=power((power(2,Rv)*sigma2/pt),-1/a);                                              %积分下限
d2=power((1/2*power(2,Rv)*sigma2/pt),-1/a);                                          %积分上限
for i=1:1:10
    pc1=1-exp(-lambdab*pi*(pt/(power(2,Rs+Rv)*sigma))^(2/a));
    pc2=1-exp(-am(i).*lambdab*pi*(pt/(power(2,Rs+Rv)*sigma))^(2/a));



    %ps1=1-exp(-lambdae*pi*((pt/(power(2,Rv)*sigma2))^(2/a)));

    k=-pi.*lambdae;
    t=power(2,Rv).*sigma2/pt;
    %fun2=2.*k.*integral(@(x)(exp(k.*power(t-power(x,-a),-2/a))-exp(k.*x.^2)).*x,d1,d2);

    %ps=1-exp(-lambdae*pi*((pt/(power(2,Rv)*sigma2))^(2/a)))+2.*k.*integral(@(x)(exp(k.*power(t-power(x,-a),-2/a))-exp(k.*x.^2)).*x,d1,d2);
    ps=1-exp(-lambdae*pi*((pt/(power(2,Rv)*sigma2))^(2/a)));
    %sum1=sum1+pm(i)*(am(i).*pc1.*ps+(1-am(i)).*pc2.*ps.*ps);
    %sum2=sum2+pm(i)*(am(i).*pc1.*(1-ps)+am(i).*(1-am(i)).*pc1.*pc1.*(1-ps).*(1-ps));
    sum2=sum2+pm(i)*(am(i)*pc1*(1-ps)+(1-am(i))*pc1*(1-ps)*am(i)*pc1*(1-ps))
end
u=sum2;

end
