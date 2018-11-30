function [ fminFit ] = fMinWeightedAvg( data,MaxTACP )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
holdingData=data;
delays=[1,3,6];
for type=1:2
    for delay=1:3
        for TACP=1:MaxTACP
  TrialsOfInterestANOVA=  holdingData.TACP==TACP & abs(holdingData.ErrWRTOptimal)<40 & holdingData{:,2}==type & holdingData.DelayTime==delays(delay) & ~isnan(holdingData.GenMean);
y=holdingData.Guess(TrialsOfInterestANOVA);
Means=holdingData.GenMean(TrialsOfInterestANOVA);
Targs=holdingData.TargAngle(TrialsOfInterestANOVA);

%Want to minimix (y real-yfit)^2  yfit=x*Means+(1-X)*Targs
yfit=@(x)(x(1)*Means+x(2)*Targs);
fun = @(x)sum((y-yfit(x)).^2);
x=fmincon(fun,[0,0],[],[],[1,1],1,[0,0],[1,1]);
fminFit(TACP,delay,type)=x(1);
        end
        
        
       
    end
end
end

