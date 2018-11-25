function xf = integrate(H,dt,x,u)
    % integration for constant H.
    xf = (expm(H(u)*dt))*x;
end