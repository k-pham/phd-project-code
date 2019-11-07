function fun = fun2DGaussian(coeffs,xymesh)
% need 5 coeffs: amplitude, xmean, xsigma, ymean, ysigma

fun = coeffs(1) * exp( - ( (xymesh(:,:,1)-coeffs(2)).^2 / (2*coeffs(3)^2) + ...
                           (xymesh(:,:,2)-coeffs(4)).^2 / (2*coeffs(5)^2) ) );

end