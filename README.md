# Tractor-Trailer-Control
This repository contains the Matlab Code used to choose and simulate a Control policy for a Tractor with a Steerable Trailer.
<ul>
  <li> main.m is the core of the program: it designs the control policy and plots the trajectory. It uses the other programs for subtasks. </li>
  <li> lie_d.m and lie_p.m implement the lie derivative of a symbolic funcion and the lie bracket between two symbolic functions, respectively. </li>
  <li> interpolate.m returns the minimum degree interpolating polynomial that fits some initial conditions [f(t_0), f'(t_0), ..., f^m(t_0)] and final conditions [f(t_f), f'(t_f), ..., f^m(t_f)].</li>
</ul>
