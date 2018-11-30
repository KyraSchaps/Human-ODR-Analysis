function [ Errfits] = LinearFitErrors( data,MaxTACP )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%separately look at linear model of errWRTOptimal as fxn of Delay, TACP


for type=1:3
    if type<3
TrialsOfInterestANOVA=  data.TACP<=MaxTACP & abs(data.ErrWRTOptimal)<40 & data{:,2}==type;
    else
        TrialsOfInterestANOVA=  data.TACP<=MaxTACP & abs(data.ErrWRTOptimal)<40;
    end
y=abs(data.ErrWRTOptimal(TrialsOfInterestANOVA));
MeanDiffs=data.MeanTargDiff(TrialsOfInterestANOVA);
TACPS=data.TACP(TrialsOfInterestANOVA);
DELAYS=data.DelayTime(TrialsOfInterestANOVA);
TYPE=data{TrialsOfInterestANOVA,2};
TABLE=table(y,MeanDiffs,TACPS,DELAYS,TYPE,'VariableNames',{'Errs','MeanDiffs','TACP','DELAYS','TYPE'});
TABLE.TYPE=categorical(TABLE.TYPE);
% TABLE.TACP=categorical(TABLE.TACP);
if type==3
p=anovan(y,{TACPS,DELAYS,TYPE},'model','interaction','varnames',{'TACPS','DELAY','TYPE'});
end
if nnz(TrialsOfInterestANOVA)>0
Errfits{type}=fitlm(TABLE,'Errs~TYPE*TACP*DELAYS');
else
    Errfits{type}=nan;
end
end

