function[Kstar,Ksat,rho,vp,vs] = gassmann(phi, S_water_0, S_water_1)

K_CLAY   = 15e+9;
K_QUARTZ = 37e+9;
F_CLAY   = 0.15;
F_QUARTZ = 1 - F_CLAY;
RHO_QUARTZ = 2650;

K_MINERAL_MATRIX = 0.5 * (F_CLAY*K_CLAY + F_QUARTZ*K_QUARTZ + ...
                          1.0 / (F_CLAY/K_CLAY + F_QUARTZ/K_QUARTZ));
                      
% Initial assumptions
Rho_initial = 2650;
Vp_initial  = 6000;
Vs_initial  = 4000;
K_initial = Rho_initial * (Vp_initial^2 - 4.0/3.0*Vs_initial^2);
G_initial = Rho_initial * Vs_initial^2;

% Fluid properties
RHO_WATER = 1000;
VP_WATER  = 1500;
K_WATER = RHO_WATER * VP_WATER^2;

RHO_OIL = 600;
VP_OIL  = 1200;
K_OIL = RHO_OIL * VP_OIL^2;

% Old saturation
K_fluid_mix = 1.0 / (S_water_0/K_WATER + (1-S_water_0)/K_OIL);

% Kstar
Kstar = (K_initial*(phi*K_MINERAL_MATRIX/K_fluid_mix + 1.0 - phi) - K_MINERAL_MATRIX) / ...
        (phi*K_MINERAL_MATRIX/K_fluid_mix + K_initial/K_MINERAL_MATRIX - 1.0 - phi);
        
% New saturation
K_fluid_mix = 1.0 / (S_water_1/K_WATER + (1-S_water_1)/K_OIL);
rho_fluid_mix = S_water_1*RHO_WATER + (1-S_water_1)*RHO_OIL;

% Ksat
Ksat = Kstar + (1.0 - Kstar/K_MINERAL_MATRIX)^2 / ...
               (phi/K_fluid_mix + (1-phi)/K_MINERAL_MATRIX - Kstar/K_MINERAL_MATRIX^2);

% New seismic parameters
grain_density = RHO_QUARTZ;
rho = grain_density*(1-phi) + rho_fluid_mix*phi;
vp = sqrt((Ksat + 4.0/3.0*G_initial) / rho);
vs = sqrt(G_initial / rho);
        




                          