%{
Breakdown of the scripts in this folder:

timevo() - A function that takes in a matrix for a differential operator
and outputs the time evolution under a small timestep for a certain number
of steps, and can accept "measured data" to acommodate data assimilation
(or generate the data).

TestScriptforabash - Supposed to be a simple model to show that timevo
works... But it doesn't yet.

Lyapunov_Code - A script that give a Lyapunov function for the Lorenz '63
System

IntroAlgorithm - A script that finds the minimal mu values for data
assimilation based on the sum-of-squares of the absolute differences.

AlgExtension - A failed script that highlights the need for a pre-defined
energy function; else, we need to solve a quadratic contraint convex
optimization problem, which is much harder than sum of squares.

abash() - Made to support timevo(), it takes in the data, operator, and
timestep and calculates the next row predicted by the model using 
second-order Adams-Bashforth.

Convg - A script that is supposed to test the data assimilation constants
and show that the model converges to the true value, given some valid
constants for mu. So far, it has not worked, and I am currently still
debugging it.

YALMOSopts - Just a list of options for YALMIP and MOSEK. Not really to be
run, just there to have the options on hand.
%}