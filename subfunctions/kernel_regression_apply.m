function yte=kernel_regression_apply(xtest,KR_model_vs)
NORMILIZE_DATA = KR_model_vs.NORMILIZE_DATA;
if NORMILIZE_DATA
    xtrain_mean      = KR_model_vs.xtrain_mean;
    xtrain_std       = KR_model_vs.xtrain_std;
    ytrain_mean      = KR_model_vs.ytrain_mean;
    ytrain_std       = KR_model_vs.ytrain_std;
    xtest_normalized = normalize_dist(xtest, xtrain_mean,xtrain_std);
    xte=xtest_normalized;
else
    xte=xtest;
end
KR_model0_ytr   = KR_model_vs.KR_model0_ytr;

Nd=size(ytrain_mean,1);
Nsample=size(xtest,2);
ydisph0=zeros(Nd,Nsample);
display('Applying Kernel Regularized Least-Squares');

h = waitbar(0,'Regression per dimension started');
for d=1:Nd
    waitbar(d/Nd,h,'Regression per dimension in progress');
    ydisph0(d,:)=IWLS_test(xte,KR_model0_ytr{d});
end
waitbar(1,h,'Regression per dimension started');
close(h);

if NORMILIZE_DATA
    ytest_normalized= ydisph0;
    yte=unnormalize_dist(ytest_normalized,ytrain_mean,ytrain_std);
else
    yte = ydisph0;
end
end
