function [o]=permutedbrunnermunzel(x,y,niter)
% only perform 2-sided test
% if niter = 0, it will perform exact test using all combinations (nCr),
% if the # of combinations > 25_Choose_11, please use random sampling,
% e.g., set niter = 100000 for 100000 iterations;
% or use the approximate test: approxbrunnermunzel.m
% ported from (and added random sampling technique)
% https://github.com/toshi-ara/brunnermunzel/blob/master/src/bm_permutation_test.f

TOLER = 10^(-10);

% columnize data
x = x(:);
y = y(:);
% remove nan if any
nonnanidx = find(~isnan(x));
x = x(nonnanidx); clear nonnanidx
nonnanidx = find(~isnan(y));
y = y(nonnanidx); clear nonnanidx

dat = [x;y];

nx = length(x); r = nx;
ny = length(y);
n = nx+ny;

if niter<1
    n_nCr = nchoosek(n,r); % if n_nCr > 4,457,400 (C(25,11)) use resampling
    fprintf(1,'\nThe number of combinations = %d\n',n_nCr);
end

[rkx, tiedadjx] = tiedrank(x);
[rky, tiedadjy] = tiedrank(y);
[rkxy, tiedadjxy] = tiedrank(dat);

mx = mean(rkxy(1:r));
my = mean(rkxy(r+1:n));
pst = (my - (n - r + 1) * 0.5) ./ r;

if niter<1
    stat_=abs(bm_permutation_stat(n, r, n_nCr, dat)); % abs: 2-sided test
    if isnan(stat_)
        return;
    end
    actualstat = stat_(1);
    %pval =length(find(stat_>=actualstat))/n_nCr;
    pval = sum(stat_>=actualstat)/n_nCr;
else
    stat_=abs(bm_resampling_stat(n, r, niter, dat)); % abs: 2-sided test
    actualstat = stat_(1);
    %pval =length(find(stat_>=actualstat))/niter;
    pval = sum(stat_>=actualstat)/niter;
end

o.pval=pval; o.actualstat=actualstat; o.allstat=stat_;
fprintf(1,'\np-value = %f',pval);
fprintf(1,'\nsample estimates of P(X<Y)+.5*P(X=Y): %f\n', pst)
return


function stat_=bm_permutation_stat(n, r, n_nCr, dat)
nx = r;
ny = n - r;
stat_ = NaN*ones(n_nCr,1);
G = zeros(n,1); % G will be used to produced logical grouping variable 

% constant values (from nx and ny) to avoid multiple calculation
const_(1) = (nx + 1) * 0.5;
const_(2) = (ny + 1) * 0.5;
const_(3) = nx * 1.0 / (nx - 1);
const_(4) = ny * 1.0 / (ny - 1);

% produce the combinations; the 1st row is actual x index
idx = nchoosek([1:n],r);
% check 1st row
if sum(idx(1,:)') ~= sum([1:r]')
    % impossible to be here :-)
    fprintf(1,'\nERROR in generating nchoosek, 1st row is wrong,')
    fprintf(1,'\nreturn to caller with stat_ has a single NaN\n');
    stat_ = NaN;
    return;
end
fprintf(1,'\nPerforming exact Brunner Munzel test in %d iterations\n',n_nCr)
% percent progress indicator, every 25%
kk = [round(n_nCr/4) round(n_nCr/2) round(3*n_nCr/4) n_nCr 0];
jj = [25 50 75 100 NaN];
ll = 1;
% start analysis (get statistics in all combinations)
for m = 1:n_nCr
    stat_(m)=calc_statistics(nx, ny, dat, const_, idx(m,:), G);
    if m == kk(ll)
        fprintf(1,'%d percent\n',jj(ll))
        ll = ll+1;
    end        
end
return

function stat_=bm_resampling_stat(n, r, niter, dat)
nx = r;
ny = n - r;
stat_ = NaN*ones(niter,1);
G = zeros(n,1); % G will be used to produced logical grouping variable 

% constant values (from nx and ny) to avoid multiple calculation
const_(1) = (nx + 1) * 0.5;
const_(2) = (ny + 1) * 0.5;
const_(3) = nx * 1.0 / (nx - 1);
const_(4) = ny * 1.0 / (ny - 1);

% the 1st stat_ uses actual x index
idx = [1:r];
stat_(1)=calc_statistics(nx, ny, dat, const_, idx, G);

fprintf(1,'\nPerforming random sampling Brunner Munzel test in %d iterations\n',niter)
% percent progress indicator, every 25%
kk = [round(niter/4) round(niter/2) round(3*niter/4) niter 0];
jj = [25 50 75 100 NaN];
ll = 1;
% start analysis (get statistics from resampled data)
for m = 2:niter
    stat_(m)=calc_statistics(nx, ny, dat, const_, randperm(n,r), G);
    if m == kk(ll)
        fprintf(1,'%d percent\n',jj(ll))
        ll = ll+1;
    end        
end
return

function stat_=calc_statistics(nx, ny, dat, const_, idx, G)
[x, y, xy] = divide_groups(nx, ny, dat, idx, G);
[rkx, tiedadjx] = tiedrank(x);
[rky, tiedadjy] = tiedrank(y);
[rkxy, tiedadjxy] = tiedrank(xy);
mx = mean(rkxy(1:nx));
my = mean(rkxy(nx+1:nx+ny));
dx = (rkxy(1:nx) - rkx - mx + const_(1)).^2;
dy = (rkxy(nx+1:nx+ny) - rky - my + const_(2)).^2;
vx = sum(dx(1:nx));
vy = sum(dy(1:ny));
v = const_(3) * vx + const_(4) * vy;
stat_ = (my - mx) ./ sqrt(v);

% start debugging
% if stat_>10
%     zz=1;
% end
% end debugging

return

function [x, y, xy] = divide_groups(nx, ny, dat, idx, G)
% G will be used to produced logical grouping variable 
% G = zeros(nx+ny,1); % moved to input to speed up.
G(idx) = 1;
Lx = logical(G);
x = dat(Lx);
y = dat(~Lx);
xy = [x;y];
return








