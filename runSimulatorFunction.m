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
minimum = 0;
times = 10;
results_mean = zeros(1, times);
results_min = zeros(1, times);

for i = 1:times
    [AvgAvail, MinAvail]= simulatorFunction(N,S,W,dlt,T,AP,pl);
    results_mean(i) = AvgAvail;
    results_min(i) = MinAvail;
end

m = mean(results_mean);
% 90% confidence interval
alfa = 0.1;
term = norminv(1-alfa/2)*sqrt(var(results_mean)/times);
fprintf('resultado = %.2e +- %.2e\n', m, term)

minimum = min(results_min);
term = norminv(1-alfa/2)*sqrt(var(results_min)/times);
fprintf('resultado = %.2e +- %.2e\n', minimum, term)

