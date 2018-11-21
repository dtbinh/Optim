import casadi.*

function xf = integrate(H,dt,x,u)
    % integration for constant H.
    xf = (expm(H(u)*dt)-1)*x
end