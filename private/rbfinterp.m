function [f,df] = rbfinterp(x, options)
tic;
phi       = options.('rbfphi');
rbfconst  = options.('RBFConstant');
nodes     = options.('x');
rbfcoeff  = (options.('rbfcoeff'))';


[dim              n] = size(nodes);
[dimPoints  nPoints] = size(x);

if (dim~=dimPoints)
  error(sprintf('x should have the same number of rows as an array used to create RBF interpolation'));
end;

f = zeros(1, nPoints);
df = zeros(dimPoints, nPoints);
% r = zeros(1, n);

for i=1:1:nPoints
    
    ds=zeros(dimPoints,1);
    dx=(x(:,i)*ones(1,n)) - nodes;
    r = sqrt(sum(dx.*dx, 1));
    
    [rbf,d_rbf]=feval(phi, r, rbfconst);
    
%     First term is constant coefficient from monomial terms
    s = rbfcoeff(n+1) + sum(rbfcoeff(1:n).*rbf);
 
	for k=1:dim
       s=s+rbfcoeff(k+n+1)*x(k,i);     % add monomial terms
%        Second term here are coefficients from monomial terms
       ds(k)=sum(rbfcoeff(1:n).*dx(k,:).*d_rbf) + rbfcoeff(k+n+1);
	end
	f(i) = s;
    df(:,i)=ds;
end

if (strcmp(options.('Stats'),'on'))
    fprintf('Interpolation at %d points was computed in %e sec\n', length(f), toc);    
end;
