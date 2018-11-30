function [ Errfits] = PropMeanLinFit( PropMeans,MaxTACP,delays)



Errfits(1:7,1)={'Intercept';'Delay'; 'TACPS'; 'Type'; 'Delay*Type'; 'Tacp*Type';'Delays*TACPS'};
y=[];
TACPs=[];
Delay=[];
Type=[];

for Condition=1:3
    y=[];
TACPs=[];
Delay=[];
Type=[];
    if Condition==3
        for type=1:2
    for delay=1:length(delays)
        for TACP=1:MaxTACP
y=[y;PropMeans(TACP,delay,type)];
TACPs=[TACPs;TACP];
Delay=[Delay;delays(delay)];
Type=[Type;type-1];
        end
    end
        end

    else
    for delay=1:length(delays)
        for TACP=1:MaxTACP
y=[y;PropMeans(TACP,delay,Condition)];
TACPs=[TACPs;TACP];
Delay=[Delay;delays(delay)];
        end

    end

        end
 
    banana=~isnan(y);
y=y(banana);
TACPs=TACPs(banana);
Delay=Delay(banana);
    if Condition==3
        Type=Type(banana);
                yfit=@(x) x(1)+x(2)*Delay+x(3)*TACPs+x(4)*Type+x(5)*Delay.*Type+x(6)*TACPs.*Type+x(7)*Delay.*TACPs;

    else
                    yfit=@(x) x(1)+x(2)*Delay+x(3)*TACPs+x(4)*Delay.*TACPs;

    end



%Want to minimix (y real-yfit)^2  yfit=x*Means+(1-X)*Targs subject to
%constraint 0<x<1
fun = @(x)sum((y-yfit(x)).^2);
if ~isnan(y);
    if Condition==3
x=fmincon(fun,[0,1,1,1,1,1,1],[],[],[],[],[0],[]);
    else
        x=fmincon(fun,[0,1,1,1],[],[],[],[],[0],[]);

    end
for i=1:length(x)
Errfits(i,Condition+1)={x(i)};
end
else
    
   for i=1:7
Errfits(i,Condition+1)={nan};
end 
end
end
end


    

