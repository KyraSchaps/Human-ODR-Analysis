function [ ErrfitsMean, ErrfitsChange] = fMinLinearfitPupilStuff( data,MaxTACP, LargeCutOff, Cutoff, SessCut, endsession,Errsize,delays )

holdingData=data;
CUTOFF = LargeCutOff;% %45; %was 75
Lcp = [[holdingData.CPMagnitude]>CUTOFF];

ErrfitsMean(1:7,1)={'Intercept';'Delay'; 'TACPS'; 'Type'; 'Delay*Type'; 'Tacp*Type';'Delays*TACPS'};
ErrfitsChange(1:7,1)={'Intercept';'Delay'; 'TACPS'; 'Type'; 'Delay*Type'; 'Tacp*Type';'Delays*TACPS'};
           
if Cutoff 
               CPCondition=Lcp;
           else
               CPCondition=ones(size(Lcp));
           end
           if endsession
               SessCondition=data.TrialInSess>(max(data.TrialInSess)-SessCut);
           else
               SessCondition=ones(size(Lcp));
        end
           
for type=1:3
    
    
    if type<3
TrialsOfInterestANOVA= SessCondition & CPCondition &  holdingData.TACP<=MaxTACP & abs(holdingData.ErrWRTOptimal)<Errsize & holdingData{:,2}==type;
        yfit=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*DELAYS.*TACPS;
    else
        TrialsOfInterestANOVA=  CPCondition &  holdingData.TACP<=MaxTACP & abs(holdingData.ErrWRTOptimal)<Errsize;
        yfit=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*TYPE+x(5)*DELAYS.*TYPE+x(6)*TACPS.*TYPE+x(7)*DELAYS.*TACPS;

    end
if nnz(TrialsOfInterestANOVA)>0
ybaseline=cell2mat(holdingData.MattMeanPupil(TrialsOfInterestANOVA));
ybaselinepast=[nan;ybaseline(1:end-1)];

ychange=cell2mat(holdingData.MattChangePupil(TrialsOfInterestANOVA));
ychangepast=[nan;ychange(1:end-1)];
ybaseline=ybaseline(~isnan(ybaseline));
ychange=ychange(~isnan(ychange));
ychangepast=ychangepast(~isnan(ychange));
ybaselinepast=ybaselinepast(~isnan(ychange));
TACPS=holdingData.TACP(TrialsOfInterestANOVA);
TACPS(TACPS==0)=1;
TACPS=TACPS(~isnan(ybaseline));
DELAYS=holdingData.DelayTime(TrialsOfInterestANOVA);
DELAYS=DELAYS(~isnan(ybaseline));
TYPE=holdingData{TrialsOfInterestANOVA,2}-1;
TYPE=TYPE(~isnan(ybaseline));


% for i=1:length(delays)
% figure(type)
% hold on
% scatter(TACPS(DELAYS==delays(i))+6*(i-1),y(DELAYS==delays(i)));
% end
% figure(type);
% boxplot(y, MaxTACP*(DELAYS)+TACPS);

% anovan(y,{TACPS,DELAYS},'model','interaction','varnames',{'TACPS','Delays'})


if type<3
            yfitbaseline=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*DELAYS.*TACPS+x(5)*ybaselinepast+x(6)*ychangepast;
            yfitchange=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*DELAYS.*TACPS+x(5)*ybaselinepast+x(6)*ychangepast+x(7)*ybaseline;

            
else
            yfit=@(x) x(1)+x(2)*DELAYS+x(3)*TACPS+x(4)*TYPE+x(5)*DELAYS.*TYPE+x(6)*TACPS.*TYPE+x(7)*DELAYS.*TACPS;

end

%Want to minimix (y real-yfit)^2  yfit=x*Means+(1-X)*Targs subject to
%constraint 0<x<1
fun1 = @(x)nansum((ybaseline-yfitbaseline(x)).^2);
fun2= @(x)nansum((ychange-yfitchange(x)).^2);
if type<3
x=fmincon(fun1,[1,1,1,1,1,1],[],[],[],[],[],[]); 
z=fmincon(fun2,[1,1,1,1,1,1,1],[],[],[],[],[],[]); 


else
% x=fmincon(fun1,[1,1,1,0,0,0,1],[],[],[],[],[],[]);
% z=fmincon(fun2,[1,1,1,0,0,0,1],[],[],[],[],[],[]);
% DELAYS=categorical(DELAYS);
% TYPE=categorical(TYPE);
% test=table(ybaseline,DELAYS,TACPS, TYPE, 'VariableNames',{'err','DELAYS','TACPS', 'TYPE'});
% fit=fitlm(test,'err~DELAYS*TYPE*TACPS')
% disp('done');
end
for i=1:length(x)
ErrfitsMean(i,type+1)={x(i)};
ErrfitsChange(i,type+1)={z(i)};
end
else
 for i=1:length(x)
ErrfitsMean(i,type+1)={nan};
ErrfitsChange(i,type+1)={nan};
end   
end    
end

    ErrfitsMean=cell2table(ErrfitsMean,'VariableNames',{'Factor','Only_Per','Only_Pre','All'});
    ErrfitsChange=cell2table(ErrfitsChange,'VariableNames',{'Factor','Only_Per','Only_Pre','All'});

end
