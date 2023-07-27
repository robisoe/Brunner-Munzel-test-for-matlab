# Brunner-Munzel-test-for-matlab
Brunner Munzel test for matlab (adapted to matlab from https://github.com/toshi-ara/brunnermunzel ). 

Files:

approxbrunnermunzel.m

permutedbrunnermunzel.m

Test data:

x=[1,2,1,1,1,1,1,1,1,1,2,4,1,1]; y=[3,3,4,3,1,2,3,1,1,5,4];

result=approxbrunnermunzel(x,y,1); % 2-sided test

 Brunner-Munzel Test Statistic = 3.137467, df = 17.682842, p-value = 0.005786
 
 sample estimates of P(X<Y)+.5*P(X=Y): 0.788961
 
 95 percent confidence interval: 0.595217, 0.982705
 

tic;result=permutedbrunnermunzel(x,y,0);toc

The number of combinations = 4457400

Performing exact Brunner Munzel test in 4457400 iterations

25 percent
50 percent
75 percent
100 percent

p-value = 0.008038

sample estimates of P(X<Y)+.5*P(X=Y): 0.788961

Elapsed time is 242.940087 seconds.



References:

Brunner, E.; Munzel, U. (2000). The nonparametric Behrens-Fisher problem: Asymptotic theory and a small-sample approximation. Biometrical Journal. 42 (1): 17–25.

Neubert, K.; Brunner, E. (2007). A studentized permutation test for the non-parametric Behrens-Fisher problem. Computational Statistics & Data Analysis. 51 (10): 5192–5204.


