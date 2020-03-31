% simulator values
N = 80;    
S = 6;      
W = 80;     
dlt = 1;   
T = 7200;  
AP = [250 1000.99965];   
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


m = sum(results_mean)/times;
minimum = min(results_min);
% 90% confidence interval
alfa = 0.1;
term_mean = norminv(1-alfa/2)*sqrt(var(results_mean)/times);
r_mean = ["result mean: ", m, " w/ confidence: ", term_mean]
term_min = norminv(1-alfa/2)*sqrt(var(results_min)/times);
r_min = ["result min: ", minimum, " w/ confidence: ", term_min]
%fprintf('resultado mean = %.2e +- %.2e\n', m, term)


%term = norminv(1-alfa/2)*sqrt(var(results_min)/times);
%fprintf('resultado min = %.2e +- %.2e\n', minimum, term)

