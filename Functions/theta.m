function out = theta(x)
global rhop10 Hs
rhop=10.^rhop10;
theta=exp(-rhop.*log(1.8.*Hs.*(0.5.*x)).^3);
theta(isinf(theta)) = 0;
out=[theta];