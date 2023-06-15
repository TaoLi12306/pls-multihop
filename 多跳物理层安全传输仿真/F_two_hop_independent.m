function u= F_two_hop_independent(Rs,Rv,Am,i,lambdab,lambdae)
%F_TWO_HOP �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
pm=zipf(10,i);%����m���������
pt=1000;     %���书��
sigma=1 ; %��������
sigma2=1;
a=4;       %˥��ָ��
Rs=Rs;      %��������
lambdab=lambdab;  %��վ�ܶ�
lambdae=lambdae;    %�����ߵ��ܶ�
Rv=Rv;        %��������
am=Am; %�洢���� 


sum1=0;
sum2=0;

d1=power((power(2,Rv)*sigma2/pt),-1/a);                                              %��������
d2=power((1/2*power(2,Rv)*sigma2/pt),-1/a);                                          %��������
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
