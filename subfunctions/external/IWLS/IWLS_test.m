function y=IWLS_test(x,model)
%
% Regression testing
% Usage:
%       y=IWLS_test(x,model)
% Input:
%    x:            d by nte test input sample matrix
%    model.alpha:  Regression parameter vector learned by IWLS
%    model.sigma:  Gaussian width chosen by cross validation
%    model.lambda: Regularization parameter chosen by cross validation
%    model.gamma:  Flattening parameter chosen by cross validation
%    model.center: Gaussian centers
% Output:
%    y:            1 by nte test output vector
%
% (c) Masashi Sugiyama, Department of Compter Science, Tokyo Institute of Technology, Japan.
%     sugi@cs.titech.ac.jp,     http://sugiyama-www.cs.titech.ac.jp/~sugi/software/uLSIF/

[d,n]=size(x);
b=size(model.center,2);
K=exp(-(repmat(sum(x.^2,1)',[1 b])+repmat(sum(model.center.^2,1),[n 1])...
        -2*x'*model.center)/(2*model.sigma^2));
y=(K*model.alpha)';
