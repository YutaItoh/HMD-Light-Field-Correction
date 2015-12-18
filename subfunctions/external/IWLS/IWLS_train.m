function model=IWLS_train(x,y,w,xte,sigma_list,lambda_list,gamma_list,b,fold)
%
% Regression training by importance-weighted least-squares
% with model selection by importance-weighted cross-validation
%
% Usage:
%       model=IWLS_train(x,y,w,xte,sigma_list,lambda_list,gamma_list,b,fold)
% Input:
%    x:              d by n training input sample matrix
%    y:              1 by n training output sample vector
%    w:              (OPTIONAL) 1 by n importance vector for training samples
%    x_te:           (OPTIONAL) d by nte test input sample matrix,
%                    which is used as Gaussian centers
%    sigma_list:     (OPTIONAL) Gaussian width
%                    If sigma_list is a vector, one of them is selected by cross validation.
%                    If sigma_list is a scalar, this value is used without cross validation.
%                    If sigma_list is empty/undefined, Gaussian width is chosen from
%                    some default canditate list by cross validation.
%    lambda_list:    (OPTIONAL) regularization parameter
%                    If lambda_list is a vector, one of them is selected by cross validation.
%                    If lambda_list is a scalar, this value is used without cross validation
%                    If lambda_list is empty, Gaussian width is chosen from
%                    some default canditate list by cross validation
%    gamma_list:    (OPTIONAL) flattening parameter
%                    If gamma_list is a vector, one of them is selected by cross validation.
%                    If gamma_list is a scalar, this value is used without cross validation
%                    If gamma_list is empty, Gaussian width is chosen from
%                    some default canditate list by cross validation
%    b:              (OPTINLAL) positive integer representing the number of kernels (default: 100)
%    fold:           (OPTINLAL) positive integer representing the number of folds
%                    in cross validation (default: 5)
%
% Output:
%      model.alpha:  Regression parameter vector learned by IWLS
%      model.sigma:  Gaussian width chosen by cross validation
%      model.lambda: Regularization parameter chosen by cross validation
%      model.gamma:  Flattening parameter chosen by cross validation
%      model.center: Gaussian centers
%
% (c) Masashi Sugiyama, Department of Compter Science, Tokyo Institute of Technology, Japan.
%     sugi@cs.titech.ac.jp,     http://sugiyama-www.cs.titech.ac.jp/~sugi/software/uLSIF/

  if nargin<2
    error('number of input arguments is not enough!!!')
  end

  [d,n]=size(x);
  [dy,ny]=size(y);
  if dy~=1
    error('y should be a horizontal vector of labels!!!')
  end
  if n~=ny
    error('sample size of x and y are different!!!')
  end
  
  if nargin<3 || isempty(w)
    w=ones(1,n);
  elseif sum(w<0)>0
    error('Importance weights must be non-negative')
  end

  if nargin<6 || isempty(lambda_list)
    lambda_list=logspace(-5,-1,9); % Candidates of regularization parameter
  elseif sum(lambda_list<0)>0
    error('regularization parameter must be non-negative')
  end

  if nargin<7 || isempty(gamma_list)
    gamma_list=linspace(0,1,11); % Candidates of flattening parameter
  elseif sum(gamma_list<0)>0
    error('flattening parameter must be non-negative')
  end
  
  if nargin<8 || isempty(b)
    b=min(100,n);
  end

  if nargin<9 || isempty(fold)
    fold=5;
  end

  %%%%%%%%%%%%%%%% Choose Gaussian kernel center `xce'
  if nargin<4
    xtmp=x;
  else
    xtmp=[x xte];
  end
  rand_index=randperm(size(xtmp,2));
  center=xtmp(:,rand_index(1:b));
  XX=repmat(sum(x.^2,1)',[1 b])+repmat(sum(center.^2,1),[n 1])-2*x'*center;

  if nargin<5 || isempty(sigma_list)
%     xscale=sqrt(median(XX(:)));
%     sigma_list=xscale*[1/10 1/5 1/2 2/3 1 1.5 2 5 10]; % Candidates of Gaussian width
    sigma_list=logspace(-2,2,9);
  elseif sum(sigma_list<=0)>0
    error('Gaussian width must be positive')
  end

 
  if length(sigma_list)==1 && length(lambda_list)==1  && length(gamma_list)==1 % need cross-validation?
    sigma_chosen=sigma_list;
    lambda_chosen=lambda_list;
    gamma_chosen=gamma_list;
  else 
    %%%%%%%%%%%%%%%% Searching Gaussian kernel width `sigma_chosen'
    %%%%%%%%%%%%%%%% and flattening parameter `lambda_chosen' 
    score_cv=zeros(length(gamma_list),length(lambda_list));
    cv_index=randperm(n);
    cv_split=floor([0:n-1]*fold./n)+1;
    
    sigma_len = length(sigma_list);
    gamma_len = length(gamma_list);
    lambda_len= length(lambda_list);
    iteration_num = sigma_len*gamma_len*lambda_len*fold;
    iteration_idx=0;
    h = waitbar(iteration_idx/iteration_num,'Cross Validation started');
    for sigma_index=1:sigma_len
      sigma=sigma_list(sigma_index);
      Ksigma=exp(-XX/(2*sigma^2));
      for gamma_index=1:gamma_len
        gamma=gamma_list(gamma_index);
        for lambda_index=1:lambda_len
          lambda=lambda_list(lambda_index);
          score_tmp=zeros(1,fold);
          for k=1:fold
            Kcvtr=Ksigma(cv_index(cv_split~=k),:);
            Kcvte=Ksigma(cv_index(cv_split==k),:);
            ycvtr=y(cv_index(cv_split~=k));
            ycvte=y(cv_index(cv_split==k));
            Kcv_w=Kcvtr.*repmat((w(cv_split~=k).^gamma)',[1 b]);
            alpha_cv=mylinsolve(Kcvtr'*Kcv_w+lambda*eye(b),Kcv_w'*ycvtr');
            score_tmp(k)=mean(((ycvte'-Kcvte*alpha_cv).^2).*w(cv_split==k)');
            % render progress bar
            iteration_idx = iteration_idx+1;
            progress=iteration_idx/iteration_num;
            waitbar(progress,h,'Cross Validation in progress...');
          end % for fold
          score_cv(gamma_index,lambda_index)=mean(score_tmp);
        end % for lambda_index
      end % for gamma_index
      [score_cv_tmp,lambda_chosen_index]=min(score_cv,[],2);
      [score_tmp(sigma_index),gamma_chosen_index]=min(score_cv_tmp);
      lambda_chosen_tmp(sigma_index)=lambda_list(lambda_chosen_index(gamma_chosen_index));
      gamma_chosen_tmp(sigma_index)=gamma_list(gamma_chosen_index);
    end % for sigma_index
    [score,sigma_chosen_index]=min(score_tmp);
    lambda_chosen=lambda_chosen_tmp(sigma_chosen_index);
    gamma_chosen=gamma_chosen_tmp(sigma_chosen_index);
    sigma_chosen=sigma_list(sigma_chosen_index);
    waitbar(1,h,'Cross Validation is done...');
    close(h);
  end %cross-validation

  K=exp(-XX/(2*sigma_chosen^2));
  K_w=K.*repmat((w.^gamma_chosen)',[1 b]);
  model.alpha=mylinsolve(K'*K_w+lambda_chosen*eye(b),K_w'*y');
  model.sigma=sigma_chosen;
  model.lambda=lambda_chosen;
  model.gamma=gamma_chosen;
  model.center=center;

