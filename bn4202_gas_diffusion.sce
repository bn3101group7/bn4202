clear;
clf;

//==============================================================================
// This functions computes the reaction term for a given species at a given
// compartment.
// \param which_species The species of gas of interest
//        1: NO
//        2: O2
// \param ind_r1 r1 (RBC | CFL)
// \param ind_r2 r2 (CFL | EC)
// \param ind_r3 r3 (EC | VW)
// \param ind_r4 r4 (VW | T)
// \param index Index pointing to the position in the domain
// \param u C_NO
// \param v P_O2
// \param lambda_core lambda_core required for calculating R_NO in RBC
// \param no_R_max R_NO_max required for R_NO and R_O2 in EC
// \param o2_Km_eNOS Constant for R_NO and R_O2 in EC
// \param lambda_vw Constant for R_NO in VW
// \param lambda_t Constant for R_NO in T
// \param Km Michaelis constant for appKm for R_O2 in VW and T
// \param no_C_ref Reference NO concentration for R_O2 in VW and T
// \param o2_Q_max_vw Q_{O2, max} for R_O2 in VW
// \param o2_Q_max_t Q_{O2, max} for R_O2 in T
// \return reaction_term The calculated reaction term
//==============================================================================
function reaction_term = GetReactionTerm(which_species, ind_r1, ind_r2,...
    ind_r3, ind_r4, index, u, v, lambda_core, no_R_max, o2_Km_eNOS,...
    lambda_vw, lambda_t, Km, no_C_ref, o2_Q_max_vw, o2_Q_max_t)
  // Start function
  select which_species
    case 1 then
      if index <= ind_r1 then
        reaction_term = lambda_core * u;
      elseif index <= ind_r2 then
        reaction_term = 0;
      elseif index <= ind_r3 then
        reaction_term = -no_R_max * v / (v + o2_Km_eNOS);
      elseif index <= ind_r4 then
        reaction_term = lambda_vw * u;
      else
        reaction_term = lambda_t * u;
      end
    case 2 then
      if index <= ind_r1 then
        reaction_term = 0;
      elseif index <= ind_r2 then
        reaction_term = 0;
      elseif index <= ind_r3 then
        reaction_term = no_R_max * v / (v + o2_Km_eNOS);
      elseif index <= ind_r4 then
        app_Km = Km * (1 + u / no_C_ref);
        reaction_term = o2_Q_max_vw * v / (v + app_Km);
      else
        app_Km = Km * (1 + u / no_C_ref);
        reaction_term = o2_Q_max_t * v / (v + app_Km);
      end
    else
      error('Wrong gas species');
  end
endfunction

//==============================================================================
// Model parameters
//==============================================================================
alpha = 1.3;  // [uM/Torr], Solubility
int_r = 25.0;  // [um], Internal radius
Km = 1;  // [Torr], Michaelis constant in the absence of NO
lambda_b = 382.5;  // [1/s], Hb scavenging at 40% Hct
lambda_t = 1;  // [1/s], Tissue scavenging (from last slide)
lambda_vw = 1;  // [1/s], Vascular wall scavenging (from last slide)
len_EC = 2.5;  // [um], Endothelial cell width
len_T = 100.0;  // [um], Tissue layer width
len_VW = 10.0;  // [um], Vessel wall width (from paper)
total_r = int_r + len_EC + len_T + len_VW;  // Total r
wss = 1.5;  // [Pa], Wall shear stress
wss_ref = 2.4;  // [Pa], Reference wall shear stress
//==============================================================================
// Gas Parameters
//==============================================================================
no_d_coeff = 3300.0;  // [um^2/s], Diffusion coefficient for NO
no_q_ref = 50.0;  // [uM/s], Reference/control NO production rate
no_C_ref = 27e-3;  // [uM], Reference NO concentration
no_R_max = wss / wss_ref * no_q_ref; // R_NO_max
o2_d_coeff = 2800.0;  // [um^2/s], Diffusion coefficient for O2
o2_Q_max_vw = 5.0;  // [uM/s], Max O2 consumption rate at vascular wall
o2_Q_max_t = 50.0;  // [uM/s], Max O2 consumption rate at tissue
o2_Km_eNOS = 4.7;  // [Torr]
o2_P = 70.0;  // [Torr], P_O2 in blood lumen
//==============================================================================
// Simulation parameters
//==============================================================================
h = 0.5;  // [um], space step
conv_tol = 1e-6;  // convergence
show_err = %T;  // Boolean to toggle display of error per iteration

// Check if space step size is appropriate
if modulo(total_r, h) > 1e-20 then, error('Wrong space step'); end

// Various compartments and domain space
nr = round(total_r / h) + 1;  // number of nodes in r direction
r = linspace(0, total_r, nr);

// Set up coefficient for R terms for clarity
r_coeff_no = h * h / no_d_coeff;
r_coeff_o2 = h * h / o2_d_coeff / alpha;

// CFL width
cfl = 3;  // [um], CFL width
lambda_core = lambda_b / 2 * (1 + int_r * int_r / (int_r - cfl) /...
    (int_r - cfl));

// Flag to check steady state condition
flag = %T;

// Various compartments and domain space
// Interface indexes
ind_r1 = (int_r - cfl) / h + 1;       // index denoting end of RBC core
ind_r2 = int_r / h + 1;               // index denoting end of vessel interior
ind_r3 = ind_r2 + len_EC / h;         // index denoting end of EC layer
ind_r4 = ind_r3 + len_VW / h;         // index denoting end of VW layer
// Number of nodes in selected compartments
nr_01 = ind_r1;                       // number of nodes for r0 < r <= r1
nr_12 = ind_r2 - ind_r1;              // number of nodes for r1 < r <= r2
nr_15 = nr - nr_01;                   // number of nodes for r1 < r <= r5

// Initialization
u = zeros(nr, 1);
u_new = u;
v = [o2_P * ones(nr_01, 1); zeros(nr_15, 1)];
v_new = v;
// Display interface posisions
printf('CFL width: %d\n', cfl);
printf('   r0    r1    r2    r3    r4    r5 \n');
printf('%5.1f %5.1f %5.1f %5.1f %5.1f %5.1f\n',...
    r(1), r(ind_r1), r(ind_r2), r(ind_r3), r(ind_r4), r($));

// Solve for an initial O2 profile
for i = ind_r1 + 1 : nr - 1
  // a term to simplify equation
  a = h / 2 / r(i);
  r_term = GetReactionTerm(2, ind_r1, ind_r2, ind_r3, ind_r4, i, u(i), v(i),...
      lambda_core, no_R_max, o2_Km_eNOS, lambda_vw, lambda_t, Km, no_C_ref,...
      o2_Q_max_vw, o2_Q_max_t);
  v(i) = 0.5 * ((1 - a) * v(i - 1) + (1 + a) * v(i + 1) - r_coeff_o2 *...
      r_term);
  // Flux BC
  v($) = v($ - 1);
end

iterations = 0;
// Start stopwatch
tic();
while flag then
  // Solve for NO
  for i = 2 : nr - 1
    // a term to simplify equation
    a = h / 2 / r(i);
    // Get NO reaction term
    r_term = GetReactionTerm(1, ind_r1, ind_r2, ind_r3, ind_r4, i, u(i),...
        v(i), lambda_core, no_R_max, o2_Km_eNOS, lambda_vw, lambda_t, Km,...
        no_C_ref, o2_Q_max_vw, o2_Q_max_t);
    u_new(i) = 0.5 * ((1 - a) * u_new(i - 1) + (1 + a) * u(i + 1) -...
        r_coeff_no * r_term);
    // Flux BC at r = 0 and r = r5
    u_new(1) = u_new(2);
    u_new($) = u_new($ - 1);
  end
  // Solve for O2
  for i = ind_r1 + 1 : nr - 1
    // a term to simplify equation
    a = h / 2 / r(i);
    // Get O2 reaction term
    r_term = GetReactionTerm(2, ind_r1, ind_r2, ind_r3, ind_r4, i, u(i),...
        v(i), lambda_core, no_R_max, o2_Km_eNOS, lambda_vw, lambda_t, Km,...
        no_C_ref, o2_Q_max_vw, o2_Q_max_t);
    v_new(i) = 0.5 * ((1 - a) * v_new(i - 1) + (1 + a) * v(i + 1) -...
        r_coeff_o2 * r_term);
    // Flux BC at r = r5
    v_new($) = v_new($ - 1);
  end
  err = max(abs(u - u_new));
  if show_err then, printf('%f\n', err); end
  if err < conv_tol then, flag = %F; end
  u = u_new;
  v = v_new;
  iterations = iterations + 1;
end
printf('Time taken: %f\nIterations: %d', toc(), iterations);
subplot(2, 1, 1);
plot(r', u);
subplot(2, 1, 2);
plot(r', v);
