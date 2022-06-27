function F = circfit(beta2)
global X Y Z timev
F =beta2(1) + beta2(2).*sqrt((beta2(3)^2 + beta2(4)^2) + X +beta2(3).*Y + beta2(4).*Z) - timev;