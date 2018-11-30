function [RandoData,Means] = SampleData(H, STD,Total )


Means=[floor(360*rand)];
RandoData=[Means+STD*randn];
 
for i=2:Total
if rand<H
    Means(i)=floor(360*rand);

else
    Means(i)=Means(i-1);
end
RandoData(i)=Means(i)+STD*randn;
end
RandoData=RandoData';
Means=Means';
end