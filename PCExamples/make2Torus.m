%R1: Outer radius
%R2: Inner radius
%N: Number of samples (will return as close as possible for an evenly
%sampled torus)
function [ X ] = make2Torus( R1, R2, N )
    ratio = R1/R2;
    k = round(sqrt(N/ratio));
    t1 = linspace(0, 2*pi, ratio*k + 1);
    t1 = t1(1:end-1);
    t2 = linspace(0, 2*pi, k + 1);
    t2 = t2(1:end-1);
    length(t1)
    length(t2)
    [theta, phi] = meshgrid(t1, t2);

    x = (R1 + R2.*cos(phi)) .* cos(theta);
    y = (R1 + R2.*cos(phi)) .* sin(theta);
    z = R2.*sin(phi);
    X = [x(:) y(:) z(:)];
end