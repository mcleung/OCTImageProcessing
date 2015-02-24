function [] = plotNormalDistr(mu, sigma, sf, xval, color)

yval = sf / sqrt(2*pi*sigma^2) * exp(-1*(xval - mu).^2 / sigma^2);

plot(xval, yval, 'color', color);