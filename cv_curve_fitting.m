function cv_curve_fitting(s1s2, cv)
% s1s2 is the s1s2 coupling interval in seconds
% cv is the conduction velocity

fun = @(X)sum(abs(cv-cv_curve(s1s2, X(1), X(2), X(3))));
X0(1) = 0.02; % starting slope, m
X0(2) = -10; % starting offset, b
X0(3) = median(s1s2); % starting breakpoint, ms

%[X,FVAL] = fminsearch(fun, X0)

lb(1) = 0.01; lb(2) = -10; lb(3) = 250;
ub(1) = 0.1;  ub(2) = 0;   ub(3) = 400;
[X, FVAL] = fmincon(fun, X0, [], [], [], [], lb, ub)

model_s1s2 = min(s1s2):max(s1s2);
model_cv = cv_curve(model_s1s2', X(1), X(2), X(3));


figure
plot(model_s1s2,model_cv, 'linewidth', 2);
hold on
plot(s1s2,cv, '.', 'markersize', 50);
legend('model', 'data', 'location', 'northeastoutside');
xlabel('S1S2 coupling interval, ms');
ylabel('Conduction velocity, m/s');
set(gca, 'xlim', [200 450], 'ylim', [0 1.5])


% get the maximum value
maxVelocity = cv_curve(X(3), X(1), X(2), X(3))


end