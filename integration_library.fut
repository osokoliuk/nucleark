-- Derive the maximum value and it's index (starts from 0) within a given array
def max_arr_idx (x_arr: []f64) : (f64, i64) =
  let (x_max, idx) =
    loop (x_max, idx) = (x_arr[0], 0)
    for i in 1..<length (x_arr) do
      if x_arr[i] > x_max then (x_arr[i], i) else (x_max, i)
  in (x_max, idx)

-- Compute an absolute value for a float, self-explanatory
def abs_value (x: f64) : f64 = if x >= 0 then x else -x

def linear_interpolation (x_arr: []f64) (y_arr: []f64) (x_eval: f64) : f64 =
  let diff = map (\x -> x - x_eval) x_arr
  let (x_max, idx) = max_arr_idx diff
  let ((x0, x1), (y0, y1)) =
    if x_max > x_eval
    then ((x_arr[idx - 1], x_max), (y_arr[idx - 1], y_arr[idx]))
    else ((x_max, x_arr[idx + 1]), (y_arr[idx], y_arr[idx + 1]))
  in y0 + (x_eval - x0) * (y1 - y0) / (x1 - x0)

--Linear interpolation formula
module type derivative = {
  type vec
  val f : {x: f64, y: f64} -> f64
  val make_ic : {x: f64, y: f64, dx: f64} -> vec
}

module type derivative_vec = {
  type vec
  type t
  type s

  val const1 : t
  val const2 : t

  val add : t -> t -> t
  val divide : t -> t -> t
  val multiply : t -> t -> t

  val n : i64
  val f [n] : {x: t, y: [n]t} -> [n]t
  val make_ic : {x: f64, y: [n]f64, dx: f64} -> vec
}

-- 4th order Runge-Kutta solver module with abstract dx/dy = f(x,y)
module runge_kutta (f_input: derivative) = {
  type t = f64
  type vec = {x: f64, y: f64, dx: f64}

  def compute_y1 {x = x0: t, y = y0: t, dx = dx: t} : (t, t, t) =
    let k1 = f_input.f {x = x0, y = y0}
    let k2 = f_input.f {x = x0 + dx / 2.0, y = y0 + dx * k1 / 2.0}
    let k3 = f_input.f {x = x0 + dx / 2.0, y = y0 + dx * k2 / 2.0}
    let k4 = f_input.f {x = x0 + dx, y = y0 + dx * k3}
    let y1 = y0 + dx / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    let x1 = x0 + dx
    in (x1, y1, dx)
}

-- 4th order Runge-Kutta solver module with abstract dx/dy = f(x,y0,y1,...,yn)
module runge_kutta_vec (f_input: derivative_vec) = {
  open f_input

  def compute_y1 [n] {x = xi: t, y = yi: [n]t, dx = dx: t} : (t, [n]t, t) =
    let k1 = f_input.f {x = xi, y = yi}
    let k2 = f_input.f {x = add xi (divide dx const1), y = map2 (\x y -> add x (multiply (divide dx const1) y)) yi k1}
    let k3 = f_input.f {x = add xi (divide dx const1), y = map2 (\x y -> add x (multiply (divide dx const1) y)) yi k2}
    let k4 = f_input.f {x = add xi dx, y = map2 (\x y -> add x (multiply dx y)) yi k3}
    let yf = map2 (\x y -> add x (multiply (divide dx const2) y)) yi (map4 (\x y z w -> add (add (add x (multiply const1 y)) (multiply const1 z)) w) k1 k2 k3 k4)
    let xf = add xi dx
    in (xf, yf, dx)
}

module type integral = {
  type t
  type s
  type vec

  val f : f64 -> f64
  val set_bounds : f64 -> f64 -> i64 -> vec
}

-- Integrate a function f(x) via a Simpsons rule with n subdivisions
-- f(x) here is an abstract function defined in another module
module simpsons (f_input: integral) = {
  type t = f64
  type s = i64

  -- Taken from https://futhark-lang.org/examples/literate-basics.html
  def linspace (n: s) (start: t) (end: t) : [n]t =
    tabulate n (\i -> start + f64.i64 i * ((end - start) / f64.i64 n))

  def compute_integral (a: t, b: t, n: s) : t =
    -- Interpolate the function over a finite grid of values between a and b
    let x_grid = linspace n a b
    let f_grid = map f_input.f x_grid
    let dx = (b - a) / f64.i64 n
    -- Sum over even and odd indices of an array via reduce
    let odd_sum = reduce (+) 0 f_grid[1::2]
    let even_sum = reduce (+) 0 f_grid[2::2]
    -- Compute the final sum
    let final_sum = 1.0 / 3.0 * dx * (f_grid[0] + 4.0 * odd_sum + 2.0 * even_sum + last (f_grid))
    in final_sum

  -- Determine the value for the number of subdivisions so that the error is sufficiently small
  -- Evaluation of the exact value for n is tricky as Futhark does not support recursive relations
  def compute_error (a: t) (b: t) (n: t) : t =
    let dx = (b - a) / n
    let grad4f_a = (vjp f_input.f a 4f64)
    -- Automatic differentiation around a
    let grad4f_b = (vjp f_input.f a 4f64)
    -- Automatic differentiation around b
    let abs_arr = map abs_value [grad4f_a, grad4f_b]
    let (max_value, _) = max_arr_idx (abs_arr)
    let err = 1.0 / 180.0 * dx ** 4 * (b - a) * max_value
    -- Upper bound for error with n subdivisions
    in err
}
