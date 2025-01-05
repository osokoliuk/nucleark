import "integration_library"
import "constants"

-- Open constants module to import all of the important quantities
open constants

-- Taken from https://futhark-lang.org/examples/literate-basics.html
def linspace (n: i64) (start: f64) (end: f64) : [n]f64 =
  tabulate n (\i -> start + f64.i64 i * ((end - start) / f64.i64 n))

-- LCDM Hubble parameter
def H (z: f64) : f64 = H0 * f64.sqrt (Omegam0 * (1 + z) ** 3 + Omegar0 * (1 + z) ** 4 + 1 - Omegam0 - Omegar0)

-- Exact form of SFRD integrand
module integrand_growth = {
  type t = f64
  type s = i64
  type vec = (f64, f64, i64)

  def set_bounds (a: f64) (b: f64) (n: i64) = (a, b, n)

  def f (x: t) : t = (Omegam0*x**(-3))**gamma_growth/x
}

module simpsons_growth = simpsons (integrand_growth)

-- Approximation for a critical overdensity treshold with growth index
-- D(z) taken from Nesseris et al. 2008
def delta_c (z: f64) : f64 =
  let bounds = integrand_growth.set_bounds 1.0 (1.0/(1.0+z)) 1000
  let Dz = f64.exp(simpsons_growth.compute_integral bounds)
  in delta_c0/Dz

def Wf (k: []f64) (R: f64) : []f64 =
  let kR = map (* R) k
  in map4 (\x y z w -> (3 * f64.sin (x) - 3 * y * f64.cos (z)) / (w) ** 3) kR kR kR kR

-- Cosmic density variance (not squared!), P(k) is obtained via CAMB python interface
def sigma_M (Mh: []f64) (Pk: []f64) : ([]f64) =
  let R = (2.0 * GN * Mh / (Omegam0 * H0 ** 2 * c ** 3)) ** (1 / 2)
  let sigma2 = map3 (\x y z -> x ** 2 / (2.0 * f64.pi ** 2) * y * z ** 2) k_arr Pk Wf (k_arr) (R) |> simpsons
  in map (\x -> x ** (1.0 / 2.0)) sigma2

-- First crossing distributions for Sheth-Tormen and Tinker halo mass functions
def first_crossing (Mh: []f64) (Pk: []f64) (kind_HMF: #ST | #Tinker) : []f64 =
  let sigma_arr = sigma Mh Pk
  let nu = map2 (\x y -> x / y) delta_c sigma_arr
  in match kind_HMF
     case #ST -> map3 (\x y z -> A_ST * (2 * x ** 2 / f64.pi) ** (1 / 2) * (1 + nu ** (-2 * p)) * f64.exp (-nu ** 2 / 2))
     case #Tinker -> map2 (\x y -> A_Tinker * ((x / b) ** (-alpha) + 1) * f64.exp (-c / y ** 2)) sigma_arr sigma_arr

-- Halo Mass Function
def HMF (Mh: []f64) (Pk: []f64) (kind_HMF: #ST | #Tinker) : []f64 =
  let log_sigma = map f64.log sigma (Mh Pk)
  let dsigma_arr = jvp log_sigma (map f64.log Mh)
  let first_crossing_arr = first_crossing Mh Pk
  in map3 (\x y z -> -rho_avg / x * y * z) Mh first_crossing_arr dsigma_arr

-- Function used in the Behroozi et al. SMHR
def f_behroozi (x: f64) (z: f64) : f64 =
  -- Generalized polynomial that describes fitting parameters, for exact values
  -- refer to the constants.fut file or Behroozi et al. (2013)
  let a = 1.0 / (1.0 + z)
  let nu = f64.exp (-4 * a ** 2)
  let polynomial {l1 = l1: f64, l2 = l2: f64, l3 = l3: f64, l4 = l4: f64} = l1 + (l2 * (a - 1.0) + l3 * z) * nu + l4 * (a - 1.0)
  let alpha = polynomial alpha_record
  let delta = polynomial delta_record
  let gamma = polynomial gamma_record
  in -f64.log10 (10 ** (alpha * x) + 1.0) + delta * (f64.log10 (1.0 + f64.exp (x))) ** gamma / (1 + f64.exp (10 ** (-x)))

-- Star Formation Efficiency for Double Power-Law, Behroozi et al. and EMERGE SMHR
def eps_star (Mh: []f64) (z: f64) (kind_SMF: #double | #behroozi | #emerge) : []f64 =
  let a = 1.0 / (1.0 + z)
  let nu = f64.exp (-4.0 * a ** 2.0)
  let polynomial {l1 = l1: f64, l2 = l2: f64, l3 = l3: f64, l4 = l4: f64} = l1 + (l2 * (a - 1.0) + l3 * z) * nu + l4 * (a - 1.0)
  let M1 = polynomial M1_record
  let eps = polynomial eps_record
  in match kind_SMF
     case #double -> map2 (\x y -> eps0 / ((x / Mh0) ** gamma_lo + (y / Mh0) ** gamma_hi)) Mh Mh
     case #behroozi -> map (\x -> 10 ** (f64.log10 (eps * M1) + f_behroozi f64.log10 (x / M1) z - f_behroozi 0.0 z))
     -- Add tabulated values
     case #emerge -> []

-- Mass Accretion Rate as per Fakhouri et al. 2010
def MAR (Mh: []f64) (z: f64) : f64 =
  map (\x -> 25.3 * (x / 1e12) ** 1.1 * (1.0 + 0.65 * z) * (H z) / H0) Mh

def SFR (Mh: f64) (z: f64) (kind_SMF: #double | #behroozi | #emerge) : f64 =
  let eps_star_arr = eps_star Mh z kind_SMF
  let MAR_arr = MAR Mh z
  in map2 (\x y -> x * Omegab0 / Omegam0 * y) eps_star_arr MAR_arr

-- Exact form of SFRD integrand
module integrand_SFRD = {
  type t = f64
  type s = i64
  type vec = (f64, f64, i64)

  def set_bounds (a: f64) (b: f64) (n: i64) = (a, b, n)

  def f (x: t) : t = linear_interpolation mass_arr integrand_arr x
}

module simpsons_SFRD = simpsons (integrand_SFRD)

-- Star Formation Rate Density
def SFRD (z: f64) (Pk: []f64) (kind_HMF: #ST | #Tinker) (kind_SMF: #double | #behroozi | #emerge) : f64 =
  -- Initialize array for halo masses with M_min ~ 10^6 [Msun], M_inf ~ 10^18 [Msun]
  let mass_arr = linspace 1000 1e6 1e18
  -- Derive dn/dM from HMF (dn/dlogM)
  let hmf_arr = map f64.log HMF (mass_arr Pk kind_HMF)
  let sfr_arr = SFR mass_arr z kind_SMF
  let integrand_arr = map2 (\x y -> x * y) hmf_arr sfr_arr
  let bounds = integrand_SFRD.set_bounds 1e6 1e18 1000
  in simpsons_SFRD.compute_integral bounds
