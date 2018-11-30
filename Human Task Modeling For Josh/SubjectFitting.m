%Fitting Data to CPPandRelia based on H rate
load('PepleParsed.mat');
delays=[2,5];


for sub=1:16
    holdingData=data.Parsed_Data{sub,1};

    for type=1:2
        for delay=1:2

  TrialsOfInterestANOVA=  holdingData{:,2}==type & holdingData.DelayTime==delays(delay) & ~isnan(holdingData.GenMean);
y=holdingData.Guess(TrialsOfInterestANOVA);
Means=holdingData.GenMean(TrialsOfInterestANOVA);
Nexts=holdingData.NextAngle(TrialsOfInterestANOVA);
Targs=holdingData.TargAngle(TrialsOfInterestANOVA);

%Want to minimize (y real-yfit)^2
Cpp1PARAM=@(x)CPPandRelia(Targs',x(1),30);

%compare responses to model predictions
fun = @(x)sum(degAngDiff(y,Cpp1PARAM(x)').^2);

%compare GM to model predictions
fun2 = @(x)sum(degAngDiff(Means,Cpp1PARAM(x)').^2);



%fmn con
x=fmincon(fun,[.5],[],[],[],[],[.05],[1]);
FitHforResponses(sub,delay,type)=x;
x=patternsearch(fun2,[.5],[],[],[],[],[.05],[1]);
FitHforGM(sub,delay,type)=x;

        end

    end

end
