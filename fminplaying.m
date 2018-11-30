%Fitting Data to CPPandRelia based on H rate
holdingData=data.Parsed_Data{7,1};
delays=[2,5];
colors={'y','m','r','g','c','b','--y','--m','--r','--g','--c','--b',':y',':m',':r',':g',':c','b:'};
j=1;
close all;
for type=1:2
    for delay=1:2

  TrialsOfInterestANOVA=  holdingData{:,2}==type & holdingData.DelayTime==delays(delay) & ~isnan(holdingData.GenMean);
y=holdingData.Guess(TrialsOfInterestANOVA);
Means=holdingData.GenMean(TrialsOfInterestANOVA);
Nexts=holdingData.NextAngle(TrialsOfInterestANOVA);
Targs=holdingData.TargAngle(TrialsOfInterestANOVA);

%Want to minimix (y real-yfit)^2  yfit=x*Means+(1-X)*Targs subject to
%constraint 0<x<1

fun = @(x)sum(degAngDiff(y,CPPandRelia(Targs',x,15)').^2);

fun3 = @(x)sum(degAngDiff(Means,CPPandRelia(Targs',x,15)').^2);
fun4 = @(x)sum(degAngDiff(Nexts,CPPandRelia(Targs',x,15)').^2);

%Make and simulate noise:
noise=15*randn(size(y));
fun2=@(x)sum(degAngDiff(CPPandRelia(Targs',.15,15)'+noise,CPPandRelia(Targs',x,15)').^2);
[Sample,MeansSamp]=SampleData(.15,15,200);
fun5=@(x)sum(degAngDiff(MeansSamp,CPPandRelia(Sample',x,30)').^2);


possibleH=[.05:.01:1];
if ~isempty(Targs)
for i=1:length(possibleH);
 errs(i)=fun(possibleH(i));   
errs2(i)=fun2(possibleH(i));
 errs3(i)=fun3(possibleH(i));   
errs4(i)=fun4(possibleH(i));
errs5(i)=fun5(possibleH(i));



end
figure(1)
hold on
ResponseErr(j)=plot(possibleH,errs,colors{(3*(type-1)+(delay))});
title('Landscape of (Response- Model)^2');
xlabel('H for model');
ylabel('Sum of Squared Err');
figure(2)
hold on
ModelErr(j)=plot(possibleH,errs2,colors{(3*(type-1)+(delay))});
title('Landscape of (Model=.15- Model)^2');
xlabel('H for model');
ylabel('Sum of Squared Err');
figure(3)
hold on
MeansErr(j)=plot(possibleH,errs3,colors{(3*(type-1)+(delay))});
title('Landscape of (Gen Mean- Model)^2');
xlabel('H for model');
ylabel('Sum of Squared Err')';
figure(4)
hold on
NextsErr(j)=plot(possibleH,errs4,colors{(3*(type-1)+(delay))});
title('Landscape of (NextTarg- Model)^2');
xlabel('H for model');
ylabel('Sum of Squared Err')';

figure(5)
hold on
NextsErr(j)=plot(possibleH,errs5,colors{(3*(type-1)+(delay))});
title('Landscape of (Means (.15=h) - Model)^2');
xlabel('H for model');
ylabel('Sum of Squared Err')';




j=j+1;



%fmn con
x=fmincon(fun,[.5],[],[],[],[],[0],[1]);
BestHFMCResponse(delay,type)=x(1);

x=fmincon(fun2,[.5],[],[],[],[],[0],[1]);
BestHFMCModel(delay,type)=x(1);
x=fmincon(fun3,[.5],[],[],[],[],[0],[1]);
BestHFMCMean(delay,type)=x(1);
x=fmincon(fun4,[.5],[],[],[],[],[0],[1]);
BestHFMCNext(delay,type)=x(1);
x=fmincon(fun5,[.5],[],[],[],[],[0],[1]);
BestHFMCSample(delay,type)=x(1);

% %pattern search 
% x=patternsearch(fun,[.5],[],[],[],[],[0],[1]);
% BestHPSResponse(delay,type)=x(1);
% 
% x=patternsearch(fun2,[.5],[],[],[],[],[0],[1]);
% BestHPSModel(delay,type)=x(1);
% x=patternsearch(fun3,[.5],[],[],[],[],[0],[1]);
% BestHPSMean(delay,type)=x(1);
% x=patternsearch(fun4,[.5],[],[],[],[],[0],[1]);
% BestHPSNext(delay,type)=x(1);


end



    
    end
end
legend([ResponseErr(1:4)],' Per2',' Per5',' Pre2', ' Pre5');
legend([ModelErr(1:4)],' Per2',' Per5',' Pre2', ' Pre5');
legend([MeansErr(1:4)],' Per2',' Per5',' Pre2', ' Pre5');
legend([NextsErr(1:4)],' Per2',' Per5',' Pre2', ' Pre5');
