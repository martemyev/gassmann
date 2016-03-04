sat = 0.2:0.1:1;
Kstar = zeros(1, size(sat, 2));
Ksat = zeros(1, size(sat, 2));
rho = zeros(1, size(sat, 2));
vp = zeros(1, size(sat, 2));
vs = zeros(1, size(sat, 2));

for i = 1:size(sat,2)
    [Kstar(i),Ksat(i),rho(i),vp(i),vs(i)] = gassmann(0.5,0,sat(i));
end

figure
plot(sat,Kstar,sat,Ksat);

figure
plot(sat,rho);

figure
plot(sat, vp);

figure
plot(sat, vs);