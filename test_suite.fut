import "integration_library"
import "constants"

open constants

-- Exact form of dx/dy = f(x,y0,y1,...,yn)
module dxdy = {
  type t = f64
  def n : i64 = 2
  type s = [n]f64
  type vec = {x: f64, y: [n]f64, dx: f64}

  def const1 : t = 2.0
  def const2 : t = 6.0

  def add (x: t) (y: t) : t = x + y
  def divide (x: t) (y: t) : t = x / y
  def multiply (x: t) (y: t) : t = x * y

  def make_ic [n] {x = x: f64, y = y: [n]f64, dx = dx: f64} : {x: f64, y: [n]f64, dx: f64} = {x, y, dx}

  def f [n] {x = x: f64, y = y: [n]f64} : [n]f64 =
    let func_arr = [x * y[0] + y[1], y[1]] :> [n]f64
    in func_arr
}

-- Convert abstract Runge-Kutta module to a particular case
module runge_kutta_over_func = runge_kutta_vec (dxdy)


-- Test vector implementation of a Runge-Kutta method
-- ==
-- entry: test_rk4
-- input { 0f64 [1f64,1f64] 1f64 }
-- output { [3.947916666666667f64, 2.708333333333333f64] }

entry test_rk4 (x0: f64) (y0: [2]f64) (dx: f64) : [2]f64 =
  let vec_ic = dxdy.make_ic {x = x0, y = y0, dx = dx}
  let (_,y_sol,_) = runge_kutta_over_func.compute_y1 vec_ic
  in y_sol
