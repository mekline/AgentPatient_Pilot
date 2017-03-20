function [h,F,p,dof,stats,BB,EE,opt]=glm(X,Y,C,M,opt,whitening)
% GLM General Linear Model estimation and hypothesis testing.
%
%   [h,F,p,dof,stat]=GLM(X,Y,C,M) estimates a linear model of the form Y = X*B + E
%   where Y is an observed matrix of response or output variables (rows are observations, columns are output variables)
%         X is an observed design or regressor matrix (rows are observations, columns are predictor variables)
%         B is a matrix of unknown regression parameters (rows are predictor variables, columns are output variables)
%         E is an matrix of unobserved multivariate normally distributed disturbances with zero mean and unknown covariance.
%   and tests a general null hypothesis of the form h = C*B*M' = 0
%   where C is matrix or vector of "predictor" contrasts (rows are contrasts, columns are predictor variables, defaults to C=eye(size(X,2)) )
% 		  M is matrix or vector of "outcome" contrasts (rows are contrasts, columns are outcome variables, defaults to M=eye(size(Y,2)) )
%
%   GLM returns the following information:
%		  h:   a matrix of estimated contrast effect sizes (h = C*B*M')
%		  F:   the test statistic value(s) (T,F,or Chi2 value, depending on whether h is a scalar, a vector, or a matrix. See below)
%		  p:   p-value of the test(s)
%		  dof: degrees of freedom
%         stat:name of statistic(s) used ('T'/'F'/'Chi2')
%
%   Additional information:
%   By default GLM will use a T, F, or a Chi2 statistic for hypothesis testing depending on the size of h=C*B*M'. The default options are (depending on the size of h):
%                  when size(h)=[1,1]      -> T statistic (note: one-sided t-test)
%	                      			          Examples of use: one-sided two-sample t-test, linear regression
%                  when size(h)=[1,Ns]     -> F statistic (note: equivalent to two-sided t-test when Ns=1)
%   							  	 		  Examples of use: Hotelling's two sample t-square test, two-sided t-test, multivariate regression
%                  when size(h)=[Nc,1]     -> F statistic (note: equivalent to two-sided t-test when Nc=1)
%					  			     		  Examples of use: ANOVA, ANCOVA, linear regression omnibus test
%                  when size(h)=[Nc,Ns]    -> Wilks' Lambda statistic, Bartlett's Chi2 approximation
%								     		  Examples of use: MANOVA, MANCOVA, multivariate regression omnibus test, likelihood ratio test
%   If the contrast matrix M is set to the string 'univariate' the options are limited to the 1st and 3rd options above (univariate T or F statistics) and the analyses are performed separately for each column of Y. 
%
%   In addition and for any combination of C or M contrast matrices the default option can be changed using the syntax GLM(X,Y,C,M,opt) where opt is one of the following character strings:
%               GLM(X,Y,C,M,'t') will perform a separate univariate T-test on each of the elements of the matrix h=C*B*M'
%  			    GLM(X,Y,C,M,'frow') will perform a separate multivariate F-test on each of the rows of the matrix h (collapsing across multiple outcome variables or outcome contrasts).
% 			    GLM(X,Y,C,M,'fcol') will perform a separate univariate F-test on each of the columns of the matrix h (collapsing across multiple predictor variables or predictor contrasts).
% 			    GLM(X,Y,C,M,'wilks') will perform a single multivariate Wilks' Lambda-test on the matrix h
%
% Examples:
%   GLM(X,Y) fits a model of the form Y=X*B and tests the hypothesis B=0
%   GLM(X,Y,C) fits a model of the form Y=X*B and tests the hypothesis C*B=0
%   GLM(X,Y,C,'univariate') fits a model of the form Y=X*B and tests the hypothesis C*B(:,i)=0 separately for each column of B (for each column of Y) 
%   GLM(X,Y,C,M) fits a model of the form Y=X*B and tests the hypothesis C*B*M'=0  
%
% Practical example:
%   % MANOVA (three groups, two outcome variables)
%   % Data preparation
%    N1=10;N2=20;N3=30;
%    Y1=randn(N1,2)+repmat([0,0],[N1,1]); % data for group 1 (N1 samples, population mean = [0,0])
%    Y2=randn(N2,2)+repmat([0,1],[N2,1]); % data for group 2 (N2 samples, population mean = [0,1])
%    Y3=randn(N3,2)+repmat([1,0],[N3,1]); % data for group 2 (N3 samples, population mean = [1,0])
%    Y=cat(1,Y1,Y2,Y3);
%    X=[ones(N1,1),zeros(N1,2); zeros(N2,1),ones(N2,1),zeros(N2,1); zeros(N3,2),ones(N3,1)];
%   % Sample data analyses
%    [h,F,p,dof,stat]=glm(X,Y,[1,-1,0;0,1,-1]); disp(['Multivariate omnibus test of non-equality of means across the three groups:']);disp([stat,'(',num2str(dof),') = ',num2str(F),'   p = ',num2str(p)]);
%    [h,F,p,dof,stat]=glm(X,Y,[1,-1,0]); disp(['Multivariate test of non-equality of means between groups 1 and 2:']);disp([stat,'(',num2str(dof(1)),',',num2str(dof(2)),') = ',num2str(F),'   p = ',num2str(p)]);
%    [h,F,p,dof,stat]=glm(X,Y,[-1,1,0],'univariate'); disp(['Univariate one-sided test of non-equality of means between groups 1 and 2 on each outcome variable:']);disp([stat,'(',num2str(dof),') = ',num2str(F(:)'),'   p = ',num2str(p(:)')]);
%

% alfnie@gmail
% 04/03

[N1,Nx]=size(X);
[N2,Ns]=size(Y);
if N1~=N2, error('wrong dimensions'); end
univariate=[];
if nargin<3 || isempty(C), C=eye(Nx);end%speye(Nx,Nx); end
if nargin<4 || isempty(M), M=speye(Ns,Ns); elseif ischar(M)&&strcmpi(M,'univariate'), M=speye(Ns,Ns); univariate=1; else Ns=rank(M); end
if nargin<5, opt=[]; end
if nargin<6, whitening=[]; end
if ~isempty(opt),switch(lower(opt)),case {'t','collapse_none','aa'},opt='AA';case {'frow','collapse_outcomes','ab'},opt='AB';case {'fcol','collapse_predictors','ba'},opt='BA';case {'chi2','wilks','collapse_all','bb'},opt='BB';otherwise,error(['Unknown option ',opt]); end; end
if isempty(univariate), univariate=(strcmp(opt,'AA')||strcmp(opt,'BA'))&&isequal(M,speye(Ns,Ns)); end

Nx=rank(X);
if isempty(whitening), dof=N1-Nx;
else X=whitening.*X; Y=whitening.*Y; whitening2=whitening.^2; dof=sum(whitening2).^2/sum(whitening2.^2)-Nx; end %Welch–Satterthwaite / correction for nonsphericity (%X=full(whitening*X);Y=full(whitening*Y);whitening2=whitening'*whitening;dof=trace(whitening2).^2/trace(whitening2.^2)-Nx;)

iX=pinv(X'*X);
B=iX*(X'*Y);
E=Y-X*B;
if univariate,EE=sparse(1:Ns,1:Ns,sum(abs(E).^2,1)); %Ne=1;
elseif size(E,2)<size(E,1), EE=M*(E'*E)*M'; else EE=E*M'; EE=EE'*EE; end %Ne=rank(full(EE)); end % "within" matrix

h=full(C*B*M');
r=C*iX*C';
if univariate,BB=sparse(1:Ns,1:Ns,sum(conj(h).*(pinv(r)*h),1)); 
else BB=h'*pinv(r)*h; end                                                % "between" matrix 
Nc0=rank(X*C');

if nargin<5 || isempty(opt), opt=[char('A'+(size(h,1)>1)),char('A'+((size(h,2)>1)&~univariate))]; else opt=upper(opt); end
switch(opt),
    case 'AA',                          % h: [1,1]  T-test
        k=                full(sqrt(diag(r)*diag(EE).'));
        F=			      real(h./max(eps,k))*sqrt(dof);     
        if exist('tcdf','file'),
        p=			      1-tcdf(F,dof);
        else
        p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Tcdf(F(idxvalid),dof); end
        end
        %dof=              dof; 
        stats=            'T';
    case 'AB',                          % h: [1,Ns] F-test
        EE=full(EE);
        F=			      real(sum((h*pinv(EE)).*conj(h),2)./diag(r))*(dof-Ns+1)/Ns;  
        if exist('fcdf','file'),
        p=			      1-fcdf(F,Ns,dof-Ns+1);
        else
        p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Fcdf(F(idvalid),Ns,dof-Ns+1); end
        end
        dof=              [Ns,dof-Ns+1];
        stats=            'F';
    case 'BA',                          % h: [Nc,1] F-test
        F=                full(real(diag(BB)./diag(EE))*dof/Nc0);      
        if exist('fcdf','file'),
        p=                1-fcdf(F,Nc0,dof);
        else
        p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Fcdf(F(idxvalid),Nc0,dof); end
        end
        dof=              [Nc0,dof];
        stats=            'F';
    case 'BB',                          % h: [Nc,Ns] Wilk's Lambda-test (Ns,dof,Nc0)
        EE=full(EE);
        F=                -(dof-1/2*(Ns-Nc0+1))*real(log(real(det(EE)./det(EE+BB))));
        if exist('chi2cdf','file'),
        p=                1-chi2cdf(F,Ns*Nc0);
        else
        p=nan+zeros(size(F));idxvalid=find(~isnan(F));if ~isempty(idxvalid), p(idxvalid)=1-spm_Xcdf(F(idxvalid),Ns*Nc0); end
        end
        dof=              Ns*Nc0; 
        stats=            'Chi2';
end
