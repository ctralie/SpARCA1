%R1: Outer radius
%R2: Inner radius
%N: Number of samples
function [ X ] = make2Torus( R1, R2, N )
	theta = 2*pi*rand(N, 1);
	phi = 2*pi*rand(N, 1);

    x = (R1 + R2.*cos(phi)) .* cos(theta);
    y = (R1 + R2.*cos(phi)) .* sin(theta);
    z = R2.*sin(phi);
    X = [x(:) y(:) z(:)];
end
