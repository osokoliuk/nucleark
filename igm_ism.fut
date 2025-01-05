import "integration_library"
import "perturbations"
import "constants"

-- Open constants module to import all of the important quantities
open constants

-- Initial Mass Function (proportional to)
def IMF (m: f64) (kind_IMF: #Kroupa | #Salpeter | #Chabrier) : f64 =
  match kind_IMF
  -- Kroupa et al. 2001 IMF (broken power-law)
  case #Kroupa ->
    if m < 0.08
    then k0 * m ** alpha0
    else if m >= 0.08 && m < 0.5
    then k1 * m ** alpha1
    else if m >= 0.5
    then k2 * m ** alpha2
    else 0
  -- Salpeter et al. 1955 IMF (single power-law)
  case #Salpeter -> m ** (-2.35)
  -- Chabrier et al. 2003 (log-normal) converted to [Mpc^-3]
  case #Chabrier ->
    if m < 1
    then A_ch * f64.exp (-(f64.log m - log center_Ch) ** 2 / (2 * sigma_Ch ** 2))
    else m ** (-1.3)

--def O_w
def O_w (z: f64) (Pk: []f64) (kind_HMF: #ST | #Tinker) (kind_SMF: #double | #behroozi | #emerge) : f64 =
-- Leftover remnant mass, taken from Iben & Tutukov 1984
  let mR: f64 = if m <= 6.8 then 0.11 * m + 0.45 else 1.5
  in 1

-- Exact form of dx/dy = f(x,y0,y1,...,yn)
module IGM_ISM_ODE = {
  type t = f64
  def n : i64 = 2
  type s = [n]f64
  type vec = {x: f64, y: [n]f64, dx: f64}

  def const1 : t = 2.0
  def const2 : t = 6.0

  def add (x: t) (y: t) : t = x + y
  def divide (x: t) (y: t) : t = x / y
  def multiply (x: t) (y: t) : t = x * y

  def make_ic [n] {x = x: f64, y = y: [n]f64, dx = dx: f64} : vec = {x, y, dx}

  def f [n] {x = x: f64, y = y: [n]f64} : [n]f64 =
    -- First and Second coupled ODE's encoded into an array
    let func_arr = [-y[0],-y[1]] :> [n]f64
    in func_arr
}
