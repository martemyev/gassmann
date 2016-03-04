phis = 0.2:0.1:1;
Kstar = zeros(1, size(phis, 2));
Ksat = zeros(1, size(phis, 2));
rho = zeros(1, size(phis, 2));
vp = zeros(1, size(phis, 2));
vs = zeros(1, size(phis, 2));

for i = 1:size(phis,2)
    [Kstar(i),Ksat(i),rho(i),vp(i),vs(i)] = gassmann(phis(i),0.1,0.5);
end

figure
plot(phis,Kstar,phis,Ksat);

figure
plot(phis,rho);

figure
plot(phis, vp);

figure
plot(phis, vs);
    