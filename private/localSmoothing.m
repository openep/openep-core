function [f_dash, df_dash] = localSmoothing(x, f, x_dash, smoothingLength, fillWith)
% ------------
% AC 8/2/21: local smoothing or piecewise interpolation, and locally
% centered moving least-squares gradient
% ------------
% x is the points at which field values f are defined. 
%   size(x) = [m, 3], where m is the number of points.
% f are the values at the points x. size(f) = [m, 1].
% x_dash are the points at which to query f. size(x_dash) = [n, 3].
% smoothingLength is a scalar smoothing length. 
%   size(smoothingLength) = [1, 1]
% f_dash is the estimated function at points x_dash. size(f_dash) =
% [n, 1]
% df_dash is the estimted function gradient at x_dash. size(df_dash) = [n,
% 3]
% ------------

f_dash = zeros(length(x_dash), 1);
df_dash = zeros(length(x_dash), 3);
[Idx, dists] = rangesearch(x, x_dash, smoothingLength);
no_smooth = cellfun(@isempty, Idx);
if strcmp(fillWith, 'nearest')
    nnIdx = knnsearch(x, x_dash(no_smooth, :));
    f_dash(no_smooth) = f(nnIdx);
else
    f_dash(no_smooth) = NaN;
end

kern = @(x, h) exp(-(x./h).^2);
% kern = @cubicSpline;

for i=1:numel(Idx)
    
    if ~isempty(Idx{i})
        
        dx = x(Idx{i},:) - x_dash(i, :);
        
        wi = kern(dists{i}', smoothingLength);
        psi = wi./sum(wi);
        
        Xi = zeros(3);
        for j = 1:size(dx, 1)
            Xi = Xi + psi(j)*dx(j, :)'*dx(j, :);
        end
        bi = pinv(Xi);
        
        f_dash(i) = sum(f(Idx{i}).*psi);
        
        for j=1:size(dx, 1)
            df_dash(i, :) = df_dash(i, :) + (f(Idx{i}(j)) - f_dash(i))*dx(j,:)*bi*psi(j);
        end
        
    end
    
end


    function w = cubicSpline(d, smoothingLength)
        q = d./smoothingLength;
        d1 = 1 - 1.5*q.^2.*(1 - q./2);
        d2 = 0.25*(2 - q).^3;
        w = zeros(size(d));
        w(q<=1) = d1(q<=1);
        w(q>1 & q<=2) = d2(q>1 & q<=2);
    end


end