function ddty = ppme(t,y)
% preditor-prey model equations
global a1 a2 a3 b1 b2 b3 c1 c2 d k
Turb = y(1);
Vzf = y(2);
Gp = y(3);
Vmf = d*Gp^2;

Q = k*t;
y1 = Turb.*Gp - a1 * Turb.^2 - a2* Vmf.^2 .* Turb - a3 * Vzf.^2 .* Turb;
y2 = b1 * Turb .* Vzf ./(1 + b2 * Vmf.^2) - b3 * Vzf;
y3 = -c1* Turb .* Gp - c2 * Gp + Q;

ddty = [y1; y2; y3];
end