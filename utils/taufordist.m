
function tau = taufordist(dist, delta, rangelen)
	dim = length(dist);
    covs = (delta * rangelen) .^ 2;
    [~, mahalaD] = mylog(zeros(1,dim), dist, covs);
    tau = 1 - chi2cdf(mahalaD, 1);
end

%% Compute the log probability density of a multivariate normal distribution
function [loglike, mahalaD] = mylog(x, means, jdiagcovs)
[n,d] = size(x);
k = size(means,1);
loglike = zeros(n, k);
mahalaD = zeros(n, k);
logDetCov = zeros(1, k);

for j = 1:k
    L = sqrt(jdiagcovs(:,:,j)); % a vector
    if  any(L < eps(max(L)) * d)
        error('Ill Conditioned Covariance.');
    end
    logDetCov(j) = sum(log(jdiagcovs(:,:,j)));
    
    Xcentered = bsxfun(@minus, x, means(j,:));
    xRinv = bsxfun(@times, Xcentered , (1 ./ L));
    
    mahalaD(:,j) = sum(xRinv.^2, 2);
    loglike(:,j) = - 0.5 * (mahalaD(:,j) + logDetCov(j) + d * log(2 * pi));
end
end
