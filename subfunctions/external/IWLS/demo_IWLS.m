% demo_IWLS.m
%
% (c) Masashi Sugiyama, Department of Compter Science, Tokyo Institute of Technology, Japan.
%     sugi@cs.titech.ac.jp,     http://sugiyama-www.cs.titech.ac.jp/~sugi/software/uLSIF/

clear all

rand('state',20);
randn('state',20);

ntr=200;
nte=1000;
check_lambda=linspace(0,1,3);
mutr=[1;1];
sigmatr=diag([0.5 0.5]);
mute=[2;2];
sigmate=diag([0.25 0.25]);
sigmay=0.1;
xtr=repmat(mutr,1,ntr)+sigmatr*randn(2,ntr);
xte=repmat(mute,1,nte)+sigmate*randn(2,nte);
ytr=sinc(xtr)+sigmay*randn(2,ntr);
yte=sinc(xte)+sigmay*randn(2,nte);

xdisp=repmat(linspace(-1,3,100),2,1);
ydisp=sinc(xdisp);
ptr_xdisp=pdf_Gaussian(xdisp,mutr,repmat(det(sigmatr),2,1));
pte_xdisp=pdf_Gaussian(xdisp,mute,repmat(det(sigmate),2,1));
w_xdisp=pte_xdisp./ptr_xdisp;

ptr_xtr=pdf_Gaussian(xtr,mutr,repmat(det(sigmatr),2,1));
pte_xtr=pdf_Gaussian(xtr,mute,repmat(det(sigmate),2,1));
w_xtr=pte_xtr./ptr_xtr;


%%%%%%%%%%%%%%%%%%%%%%%%% Estimating density ratio
[wh_xtr,wh_xdisp,uLSIF_sigma_chosen,uLSIF_lambda_chosen]=uLSIF(xtr,xte,xdisp);
%[wh_xtr]=uLSIF(xtr,xte);
disp(sprintf('uLSIF: sigma = %g, lambda = %g',uLSIF_sigma_chosen,uLSIF_lambda_chosen))

%%%%%%%%%%%%%%%%%%%%%%%%% Adaptive Importance-Weighted Kernel Regularized Least-Squares
[IWLS_model]=IWLS_train(xtr,ytr,wh_xtr,xte);
disp(sprintf('IWLS: sigma = %g, lambda = %g, gamma = %g'...
             ,IWLS_model.sigma,IWLS_model.lambda,IWLS_model.gamma))
%[yteh]=IWLS_test(xte,model);
[ydisph]=IWLS_test(xdisp,IWLS_model);

%%%%%%%%%%%%%%%%%%%%%%%%% Plain Kernel Regularized Least-Squares
[IWLS_model0]=IWLS_train(xtr,ytr,[],xte,[],[],0);
disp(sprintf('LS: sigma = %g, lambda = %g'...
             ,IWLS_model0.sigma,IWLS_model0.lambda))
[ydisph0]=IWLS_test(xdisp,IWLS_model0);


figure(1);clf;hold on
set(gca,'FontName','Helvetica')
set(gca,'FontSize',12)
plot(xdisp,ptr_xdisp,'b-','LineWidth',2)
plot(xdisp,pte_xdisp,'k-','LineWidth',2)
plot(xdisp,w_xdisp,'r-','LineWidth',2)
plot(xdisp,wh_xdisp,'g-','LineWidth',2,'Color',[0 200 0]/255)
plot(xtr,wh_xtr,'bo','LineWidth',1,'MarkerSize',8)
legend('p_{tr}(x)','p_{te}(x)','w(x)','w-hat(x)','w-hat(x^{tr})',2)
xlabel('x')
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 12 9]);
print('-dpng','density')


figure(2);clf;hold on
set(gca,'FontName','Helvetica')
set(gca,'FontSize',12)
plot(xdisp,ydisp,'r-','LineWidth',2)
plot(xdisp,ydisph,'-','LineWidth',2,'Color',[0 200 0]/255)
plot(xdisp,ydisph0,'m-','LineWidth',2)
plot(xtr,ytr,'bo','LineWidth',1)
plot(xte(1:min(100,nte)),yte(1:min(100,nte)),'kx','LineWidth',1)
legend('f(x)','f-hat_{IWLS}(x)','f-hat_{LS}(x)','Training','Test',1)
xlabel('x')
set(gcf,'PaperUnits','centimeters');
set(gcf,'PaperPosition',[0 0 12 9]);
print('-dpng','function')

