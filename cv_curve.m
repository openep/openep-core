function y = cv_curve(s1s2, m, b, c)

% m - is the gradient of the slope part of the curve
% b - is the intersect of the y part of the curve
% c - is the point on the x axis at which the highest part fo the curve is
% reached
% x - is the range


y = m * s1s2 + b;
y(s1s2>c) = m * c + b;

end