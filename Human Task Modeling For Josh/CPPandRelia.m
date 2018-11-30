function [B] = CPPandRelia( Data, Hrate, STD )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%Data should be 1xnDatapoints array
%From Rational regulation of learning dynamics by pupil-linked arousal
%systems" I think
%Reduced Baysian Delta Rule Model
%%
B=nan(1, length(Data));
B(1)=Data(1,1);
PCha=nan(1, length(Data));
Relia=nan(1, length(Data));
Relia(1)=.5;
alpha=nan(1,length(Data));
delta=nan(1,length(Data));
sigTot=nan(1,length(Data));
H=Hrate;

Modifier=[0,360,-360];

i=1;
 sigTot(i)=STD^2+Relia(i)*STD^2/(1-Relia(i));  %%Prob with normpdf b/c not circular **** Equation 6
    [PDrawn,index]=max([normpdf(Data(1,i),B(1,i),sigTot(1,i)),normpdf(Data(1,i)+360,B(1,i),sigTot(1,i)),normpdf(Data(1,i)-360,B(1,i),sigTot(1,i))]);
    PCha(i)=(1/360*Hrate)/(1/360*Hrate+PDrawn*(1-H));  %**** equation 5
    DataCorr=Data(1,i)+Modifier(index);   %Modifier to acount for circularity
    ThirdTerm(i)=(DataCorr*Relia(i)+B(i)*(1-Relia(i))-DataCorr)^2 ; % **** Part of equation 7, specifically the third multiplier of the third part in the numerator
    Relia(i+1)=(STD^2*PCha(i)+(1-PCha(i))*(Relia(i)*STD^2)+PCha(i)*(1-PCha(i))*ThirdTerm(i))/(STD^2*PCha(i)+(1-PCha(i))*(Relia(i)*STD^2)+PCha(i)*(1-PCha(i))*ThirdTerm(i)+STD^2);  %***All of equation 7
    delta(i)=-rad2deg(angdiff(deg2rad(Data(1,i)),deg2rad(B(1,i)))); % equation 3
    alpha(i)=Relia(i)+(1-Relia(i))*PCha(i); %equation 4

for i=2:length(Data)
    sigTot(i)=STD^2+Relia(i)*STD^2/(1-Relia(i));  %%Prob with normpdf b/c not circular **** Equation 6
    [PDrawn,index]=max([normpdf(Data(1,i),B(1,i-1),sigTot(1,i)),normpdf(Data(1,i)+360,B(1,i),sigTot(1,i)),normpdf(Data(1,i)-360,B(1,i),sigTot(1,i))]);
    PCha(i)=(1/360*Hrate)/(1/360*Hrate+PDrawn*(1-H));  %**** equation 5
    DataCorr=Data(1,i)+Modifier(index);   %Modifier to acount for circularity
    ThirdTerm(i)=(DataCorr*Relia(i)+B(i-1)*(1-Relia(i))-DataCorr)^2 ; % **** Part of equation 7, specifically the third multiplier of the third part in the numerator
    Relia(i+1)=(STD^2*PCha(i)+(1-PCha(i))*(Relia(i)*STD^2)+PCha(i)*(1-PCha(i))*ThirdTerm(i))/(STD^2*PCha(i)+(1-PCha(i))*(Relia(i)*STD^2)+PCha(i)*(1-PCha(i))*ThirdTerm(i)+STD^2);  %***All of equation 7
    delta(i)=-rad2deg(angdiff(deg2rad(Data(1,i)),deg2rad(B(1,i-1)))); % equation 3
    alpha(i)=Relia(i)+(1-Relia(i))*PCha(i); %equation 4
    B(i)=B(i-1)+alpha(i)*delta(i);  %equation 3
     B(B>360)=mod(B(B>360),360);
     B(B<=0)=B(B<=0)+360;
    
end
% i=i+1;
% sigTot(i)=STD+Relia(i)*STD/(1-Relia(i));
% PDrawn=max([normpdf(Data(1,i),B(1,i),sigTot(1,i)),normpdf(Data(1,i)+360,B(1,i),sigTot(1,i)),normpdf(Data(1,i)-360,B(1,i),sigTot(1,i))]);
% 
% PCha(i)=(1/360*Hrate)/(1/360*Hrate+PDrawn*(1-H));
% 
% EndBelief=B(i);
% EndPCha=PCha(i);
% EndRelia=Relia(i);

% Data(4,:)=B(:);
% Data(5,:)=PCha(:);
% Data(6,:)=Relia(:);
% Data(7,:)=sigTot(:);

end

