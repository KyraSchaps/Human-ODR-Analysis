function [ GM, Samples,CP ] = GenFakeCPData2( num, H, STD )
%UNTITLED6 Summary of this function goes here
%   This function genrates #=num of data points according to the H rate
%   process with Hrate=H and STD=STD.  

%   Return GM=Generative Mean, Samples from process, and the CP=TACP

%Set up the matrices and pick initial GM, Sample, and TACP
GM=[];
Samples=[];
GM(1)=ceil(360*rand); 
Samples(1)=STD*randn+GM(1);
CP(1)=1;

for i=2:num
 
%Check to see if its a Change point
if rand<H
    GM(i)=ceil(360*rand);
    CP(i)=1;
else 
    GM(i)=GM(i-1);
    CP(i)=CP(i-1)+1;
end

%pick sample from current GM
Samples(i)=STD*randn+GM(i);


end
Samples=Samples';
GM=GM';
end

