%      Matlab program for “Gassmann fluid substitution”
%      Author: Dhananjay Kumar (September 8, 2005)
%      Chevron Energy Technology Company, California, USA
%      References: Wang (2001), Batzle and Wang (1992), Geophysics
%      AIM: Model fluid properties for brine sand, oil sand and gas sand
%      Desired properties: P- and S-wave velocities, and density
%      Assumptions: mineral is a mixture of quartz and Clay
%      And as in the Gassmann theory (e.g., homogeneous fluid, isotropy)
%      input:  rho_o = reference density of oil ( 42 deg API)
%      GOR = gas-to-oil ratio (160 l/l)
%      rho_g = specific gravity of gas (0.9 API)
%      T = Temperature (150 deg C)
%      P = Pressure (3200 psi)
%      S = salinity (3800 ppm)
%      phi = porosity (0.20)
%      VSH = volume shale log (0.20)
%      isw = SWT: initial water saturation from log (0.40)
%       tsw = target water saturation (1.00)
%      ifluid = type of initial hydrocarbon (Gas, Oil)
%      fluid  = type of output fluid (Brine, Gas, Oil)
%      vp = P-wave velocity from log (ft/s) - insitu / original
%      vs = S-wave velocity from log (ft/s) - insitu / original
%      rho = Bulk density from log (gm/cc)  - insitu / original
%      Output: vp_sat = P-wave velocity after fluid subs (ft/s)
%      vs_sat = S-wave velocity after fluid subs (ft/s)
%      rho_sat = Density after fluid substitution (gm/cc)
%      How to run: check all input and enter file name on the matlab prompt
%      NOTE: if hydroc. is oil, it contains some dissolved gas (defined by GOR)
%        if desired fluid is oil or gas, it contains water (defined by tws)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters (use defined)
%
rho_o = 42;             % Oil gravity (deg API)
GOR = 160.0;            % GOR (L/L)
rho_g = 0.9;            % Gas gravity (API)
T = 150.00;             % Temperature (0 C)
P = 3200.00;            % Pressure (psi)
S = 3800;               % water salinity (ppm)
phi = 0.20;             % porosity (in fraction)
vsh = 0.20;             % Vsh (volume shale in fraction)
isw = 0.40;             % initial water saturation (SW)
tsw = 1.00;             % target water saturation (in fraction)
ifluid = 1;             % initial hydrocarbon is 1(oil), 2(gas)
fluid = 1;              % Desired fluid is 1(brine), 2(oil) 3(gas)
vp = 11000.0;           % ft/s  - from log (initial value)
vs = 6500.0;            % ft/s  - from log (initial value)
rho = 2.2;              % gm/c -  from log (initial value)
%
% Fixed parameters (e.g., Mavko et al., 1998)
%
k_clay = 20.9;
% Bulk mod (GPa)
k_qtz = 36.6;
rho_clay = 2.58;        % gm/cc
rho_qtz = 2.65;
%
% some applied properties
%
div_mill = 1/1000000;   % factor used to divide by million
fs2kms = 0.000305;      % factor for ft/s to km/s conversion
kms2fs = 3280.84;       % factor for km/s to ft/s conversion
v_clay = vsh*0.70;      % Assumption: V_clay = 70% ofVSH
v_qtz= 1-v_clay;        % quartz fraction in mineral
ish = 1-isw;            % initial hydrocarbon saturation
tsh = 1-tsw;            % final hydrocarbon saturation
rho_o = 141.5/(rho_o+131.5);  % oil gravity in gm/cc (from API)
P = P*6.894757*0.001;   % Press in MPa (from Psi)
S = S*div_mill;         % salinity as weight fraction
vp = vp*fs2kms;         % ft/s to km/s
vs = vs*fs2kms;         % ft/s to km/s
%
% Step 1: Matrix properties (using VRH averaging, equation 6)
%
k_voigt = v_clay*k_clay + v_qtz*k_qtz;
k_reuss = 1/(v_clay/k_clay + v_qtz/k_qtz);
k_matrix = 0.5*(k_voigt + k_reuss);            % GPa
rho_matrix = v_clay*rho_clay+v_qtz*rho_qtz;    % gm/cc
%
% Step 2: water/brine properties (Equations 10 and 11)
%
w(1,1) = 1402.85;         w(1,3) = 3.437*10^(-3);    % Table 2
w(2,1) = 4.871;           w(2,3) = 1.739*10^(-4);
w(3,1) = -0.04783;        w(3,3) = -2.135*10^(-6);
w(4,1) = 1.487*10^(-4);   w(4,3) = -1.455*10^(-8);
w(5,1) = -2.197*10^(-7);  w(5,3) = 5.230*10^(-11);
w(1,2) = 1.524;           w(1,4) = -1.197*10^(-5);
w(2,2) = -0.0111;         w(2,4) = -1.628*10^(-6);
w(3,2) = 2.747*10^(-4);   w(3,4) = 1.237*10^(-8);
w(4,2) = -6.503*10^(-7);  w(4,4) = 1.327*10^(-10);
w(5,2) = 7.987*10^(-10);  w(5,4) = -4.614*10^(-13);
sum = 0;
for i=1:5
    for j=1:4
        sum = sum+w(i,j)*T^(i-1)*P^(j-1);
    end
end
v_water = sum;
v1 = 1170-9.6*T+0.055*T*T-8.5*10^(-5)*T*T*T+2.6*P-0.0029*T*P-0.0476*P*P;
v_brine = v_water+S*v1+S^1.5*(780-10*P+0.16*P*P)-1820*S*S; % m/s
r1 = 489*P-2*T*P+0.016*T*T*P-1.3*10^(-5)*T*T*T*P-0.333*P*P-0.002*T*P*P;
rho_water=1+10^(-6)*(-80*T-3.3*T*T+0.00175*T*T*T+r1);
r2 = 300*P-2400*P*S+T*(80+3*T-3300*S-13*P+47*P*S);
rho_brine = rho_water+0.668*S+0.44*S*S+10^(-6)*S*r2;   % gm/cc (held const)
k_brine = rho_brine*v_brine*v_brine*div_mill;          % GPa (held const)
%
% Step 3: Initial Hydrocarbon properties (Equations 32 to 35)
%
if ifluid == 1      %’Oil’ Oil by default contains gas also
    B0 = 0.972+0.00038*(2.495*GOR*sqrt(rho_g/rho_o)+T+17.8)^1.175;
    rho_ps = rho_o/((1+0.001*GOR)*B0);
    rho_s = (rho_o+0.0012*GOR*rho_g)/B0;
    r1 = rho_s+(0.00277*P-1.71*0.0000001*P*P*P)*(rho_s-1.15)^2+3.49*0.0001*P;
    rho_hyc = r1/(0.972+3.81*0.0001*(T+17.78)^1.175); % gm/cc (will change)
    v = 2096*sqrt(rho_ps/(2.6-rho_ps))-3.7*T+4.64*P+0.0115*(sqrt(18.33/rho_ps-16.97)-1)*T*P;
    k_hyc = rho_hyc*v*v*div_mill;                     % GPa (will change)
else                    %’gas’ : means no OIL only gas is present
    R = 8.314;         % gas constant (eqn, same as in step7 for fluid == 3)
    Ta = T+273.15;
    Ppr = P/(4.892-0.4048*rho_g);
    Tpr = Ta/(94.72+170.75*rho_g);
    E1 = exp(-Ppr^1.2/Tpr*(0.45+8*(0.56-1/Tpr)^2));
    E = 0.109*(3.85-Tpr)^2*E1;
    Z1 = 0.03+0.00527*(3.5-Tpr)^3;
    Z = Z1*Ppr+0.642*Tpr-0.007*Tpr^4-0.52+E;
    rho_hyc = 28.8*rho_g*P/(Z*R*Ta);
    dz_dp = Z1+0.109*(3.85-Tpr)^2*E1*(-1.2*Ppr^0.2/Tpr*(0.45+8*(0.56-1/Tpr)^2));
    yo = 0.85+5.6/(Ppr+2)+27.1/(Ppr+3.5)^2-8.7*exp(-0.65*(Ppr+1));
    k_hyc = P*yo/1000*1.0/(1-Ppr/Z*dz_dp);    % GPa
end
%
% Step 4: Fluid properties(initial insitu model, equations 30and 31)
%
k_fl = 1/(isw/k_brine+ish/k_hyc);
rho_fl = isw*rho_brine+ish*rho_hyc;
%
% Step 5: Insitu original moduli (for saturated – insitu rock, equations 4 and 5)
%
dens_poros = 0; % 1 (use porosity to est initial density), 0 (use input log)
if dens_poros == 1
    rho = phi*rho_fl + (1-phi)*rho_matrix;
end
k_sat = rho*(vp*vp-vs*vs*4/3);   % GPa (will change in step 9)
g = rho*vs*vs;  % GPa (held constant)
%
% Step 6: Porous frame properties (rewrite Gassmann eqn, equation 36)
%
k1 = k_sat*(phi*k_matrix/k_fl+1-phi)-k_matrix;
k2 = phi*k_matrix/k_fl+k_sat/k_matrix-1-phi;
k_frame = k1/k2;                             % GPa (held constant)
%
% Step 7: select the type of output fluid, cal hyc/fluid prop (equations 32 to 35)
%
if fluid == 1         %’Brine’
%   disp(‘Changing fluid to brine’)
elseif fluid == 2   %’Oil’
%    disp(‘Changing fluid to Oil [with dissolved gas] with TWS brine’)
    B0 = 0.972+0.00038*(2.495*GOR*sqrt(rho_g/rho_o)+T+17.8)^1.175;
    rho_ps = rho_o/((1+0.001*GOR)*B0);
    rho_s = (rho_o+0.0012*GOR*rho_g)/B0;
    r1 = rho_s+(0.00277*P-1.71*0.0000001*P*P*P)*(rho_s-1.15)^2+3.49*0.0001*P;
    rho_hyc = r1/(0.972+3.81*0.0001*(T+17.78)^1.175);    % gm/cc (will change)
    v = 2096*sqrt(rho_ps/(2.6-rho_ps))-3.7*T+4.64*P+0.0115*(sqrt(18.33/rho_ps-16.97)-1)*T*P;
    k_hyc = rho_hyc*v*v*div_mill;                        % GPa (will change)
elseif fluid == 3   %’Gas’
%    disp(‘Changing fluid to Gas with TWS brine’)
    R = 8.314;           % gas constant
    Ta = T+273.15;
    Ppr = P/(4.892-0.4048*rho_g);
    Tpr = Ta/(94.72+170.75*rho_g);
    E1 = exp(-Ppr^1.2/Tpr*(0.45+8*(0.56-1/Tpr)^2));
    E = 0.109*(3.85-Tpr)^2*E1;
    Z1 = 0.03+0.00527*(3.5-Tpr)^3;
    Z = Z1*Ppr+0.642*Tpr-0.007*Tpr^4-0.52+E;
    rho_hyc = 28.8*rho_g*P/(Z*R*Ta);
    dz_dp=Z1+0.109*(3.85-Tpr)^2*E1*(-1.2*Ppr^0.2/Tpr*(0.45+8*(0.56-1/Tpr)^2));
    yo = 0.85+5.6/(Ppr+2)+27.1/(Ppr+3.5)^2-8.7*exp(-0.65*(Ppr+1));
    k_hyc = P*yo/1000*1.0/(1-Ppr/Z*dz_dp);    % GPa
end
%
% Step 8: Fluid properties (target saturation) and saturated rock density (equations 30 and 31)
%
k_fl = 1/(tsw/k_brine + tsh/k_hyc);
rho_fl = tsw*rho_brine + tsh*rho_hyc;
rho_sat = phi*rho_fl+(1-phi)*rho_matrix                 % gm/cc (OUTPUT)
%
% Step 9: Gassmann Saturated bulk modulus (equation 3)
%
k1 = phi/k_fl+(1-phi)/k_matrix-k_frame/(k_matrix*k_matrix);
k_sat_new = k_frame + ((1-k_frame/k_matrix)^2)/k1;
%
% Step 10: Seismic velocity after fluid substitution (equations 1 and 2)
%
vp_sat = sqrt((k_sat_new+g*4/3)/rho_sat)*kms2fs         % ft/s (OUTPUT)
vs_sat = sqrt(g/rho_sat)*kms2fs                         % ft/s (OUTPUT)
%%%%%%%%%%%%% end of matlab code %%%%%%%%%%%%%%%%%%