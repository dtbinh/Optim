function xf = euler( ode, h, x, u )
  xf = x + h*ode(x,u);
end

