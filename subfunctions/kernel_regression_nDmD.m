function KR_model_vs=kernel_regression_nDmD(xtrain,ytrain,misc)
if nargin<3
    misc='';
end
NORMILIZE_DATA=true;
%NORMILIZE_DATA=false;
USE_IWLS =false;
REVESE =1;
if REVESE
    tmp = xtrain;
    xtrain = ytrain;
    ytrain = tmp;
end

% normalize input data
% xtest = UVST_vs_test;
% xtest_normalized=normalize(xtest, xtrain_mean,xtrain_std);

display('Regression by Kernel Regularized Least-Squares');
if NORMILIZE_DATA
    xtrain_mean=mean(xtrain,2);
    xtrain_std=std(xtrain')';
    xtrain_normalized = normalize_dist(xtrain, xtrain_mean,xtrain_std);
    ytrain_mean=mean(ytrain,2);
    ytrain_std=std(ytrain')';
    ytrain_normalized = normalize_dist(ytrain, ytrain_mean,ytrain_std);
    xtr=xtrain_normalized;
    ytr=ytrain_normalized;
else
    xtr=xtrain;
    ytr=ytrain;
end
Nd=size(ytr,1);
KR_model0_ytr=cell(Nd,1);
if USE_IWLS == false
    %%%%%%%%%%%%%%%%%%%%%%%%% Plain Kernel Regularized Least-Squares
    x=xtr;
    lambda_list=logspace(-8,-1,3); % Candidates of regularization parameter% ISAMR 2015
    lambda_list=logspace(-9,-1,6); % Candidates of regularization parameter
    %lambda_list=1e-6;
    
    n=size(x,2);
    b=min(200,n);% 150;
    %xtmp=[xtr xte];
    xtmp=xtr;
    rand_index=randperm(size(xtmp,2));
    center=xtmp(:,rand_index(1:b));
    XX=repmat(sum(x.^2,1)',[1 b])+repmat(sum(center.^2,1),[n 1])-2*x'*center;
    xscale=sqrt(median(XX(:)));
    %
    if strcmp(misc,'static_xscale')
        xscale =50;%10; %50 for the optical distortion
    end
     %xscale =0.8;% This seemed to be 50 for ISMAR2015, 20150531
    sigma_list=xscale*[1/10 1/5 1/2 2/3 1 1.5 2 5 10]; % Candidates of Gaussian widt
%sigma_list=10
    for d=1:Nd
        KR_model0_ytr{d}=IWLS_train(xtr,ytr(d,:),[],[],sigma_list,lambda_list,0,b);
    end
    KR_model_vs.lambda_list=lambda_list;
    KR_model_vs.sigma_list=sigma_list;
else
    %%%%%%%%%%%%%%%%%%%%%%%%% Adaptive Importance-Weighted Kernel Regularized Least-Squares
    [wh_xtr]=uLSIF(xtr,xte);
    for d=1:Nd
        KR_model0_ytr{d}=IWLS_train(xtr,ytr(d,:),[],[],sigma_list,lambda_list,0,base);
    end
end
for d=1:Nd
    fprintf('LS_ytr%d: sigma = %g, lambda = %g\n',d,KR_model0_ytr{d}.sigma,KR_model0_ytr{d}.lambda)
end
% [ydisph0_ytr1]=IWLS_test(xte,KR_model0_ytr1);
% [ydisph0_ytr2]=IWLS_test(xte,KR_model0_ytr2);
% ydisph0=[ydisph0_ytr1;ydisph0_ytr2];
% if NORMILIZE_DATA
%     ytest_normalized= ydisph0;
%     yte=unnormalize(ytest_normalized,ytrain_mean,ytrain_std);
% else
%     yte= ydisph0;
% end

KR_model_vs.xtr=xtr;
KR_model_vs.ytr=ytr;
if NORMILIZE_DATA
    KR_model_vs.ytrain_mean=ytrain_mean;
    KR_model_vs.ytrain_std=ytrain_std;
    KR_model_vs.xtrain_mean=xtrain_mean;
    KR_model_vs.xtrain_std=xtrain_std;
end
KR_model_vs.USE_IWLS = USE_IWLS;
KR_model_vs.NORMILIZE_DATA=NORMILIZE_DATA;
KR_model_vs.KR_model0_ytr=KR_model0_ytr;
end
