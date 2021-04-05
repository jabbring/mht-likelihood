function [probs, grad] = numinvlap(ltfun,par,y,varargin)
% ///////////////////////////////////////////////////////////////////////
% Numerical inversion of Laplace transform of distribution MHT model
% with model-specific laplace transform function 'ltfun'
%
% Author: Tim Salimans (minor edit by Jaap Abbring)
%
% INPUT:
% - ltfun, function handle to calculate exp(psi(s)*t) times laplace transform
% - par, k x 1 parameter vector, first element is log variance of BM component
% - y, n x 1 vector of durations
% - varargin, number of arguments to be passed to 'ltfun'
%
% OUTPUT:
% - probs, n x 1 vector of hitting time densities
% - grad, n x k jacobian of 'probs'
%
% ///////////////////////////////////////////////////////////////////////

% parameters for numerical inversion
A=22;
l=1;
N=9;
M=25;

% % as in Rogers(2000)
% A=22;
% l=1;
% N=6;
% M=15;

% pre-calculate
T=N+M;
DT=2*T+1;

% retrieve variance
var=exp(par(1));

% evaluate BM laplace exponent and its derivative
z=(1./y)*((pi*1i/l)*(-T:T)+0.5*A/l);
zvar=1+2*var*z;
sqrz=sqrt(zvar);
lap0inv=(sqrz-1)/var;
lap0invgrad=-1./sqrz;

% calculate exp(psi(Lambda_BM(z))*t+I(~cens).log(t)) x Laplace transform
if nargout==1
    [lapprod,lapgrad]=ltfun(par,lap0inv,y,varargin{:});
else
    % and derivatives
    [lapprod,lapgrad,Glapprod,Glapgrad]=ltfun(par,lap0inv,y,varargin{:});
end

% matrix of integration terms
pmat=repmat((1/(2*l))./y,1,DT);
S=pmat.*real(lapgrad.*lap0invgrad.*lapprod);

% derivatives of integration terms, if asked
if nargout==2
    G=cell(length(par)+1,1);
    for i=1:length(G)
        G{i}=zeros(size(pmat));
        if ~isempty(Glapgrad{i})
            G{i}=G{i}+pmat.*lap0invgrad.*Glapgrad{i}.*lapprod;
        end
        if ~isempty(Glapprod{i})
            G{i}=G{i}+pmat.*lap0invgrad.*lapgrad.*Glapprod{i};
        end
    end
    
    % gradient of lap0invgrad to log variance
    G{1}=G{1}+pmat.*(z./zvar.^1.5).*lapgrad.*lapprod*var;
    
    % gradient of lap0inv to log variance
    G{1}=G{1}+G{end}.*(z./sqrz-lap0inv);
    
end

% weights to be used in Euler summation
w=zeros(DT,1);
nmrs=abs(-T:T);
for k=0:M
    ind=(nmrs<=N+k);
    w(ind)=w(ind)-nchoosek(M,k)/(2^M);
end

% calculate probabilities
probs=S*w;

% calculate gradient
if nargout==2
    grad=zeros(length(y),length(par));
    for i=1:length(par)
        grad(:,i)=real(G{i})*w;
    end
end