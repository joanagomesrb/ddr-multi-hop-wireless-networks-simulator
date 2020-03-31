% simulator values
N = 60;    
S = 3;      
W = 60;     
dlt = 1;   
T = 7200;  
AP = [250 100];   
pl = 0;  

% mean and minimum values 
m = 0;
min = 0;
times = 10;

for i = 1:times
    [AvgAvail, MinAvail]= simulatorFunction(N,S,W,dlt,T,AP,pl);
    m = m + AvgAvail;
    min = min + MinAvail;
end

mean = (m/times)*100
minimum = (min/times)*100

