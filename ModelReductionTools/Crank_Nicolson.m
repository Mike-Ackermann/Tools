function [A,B] = Crank_Nicolson(A,B,delta_t)
% For converting a continuous time system to discrete time via the
% Crank-Nicholson method.

% A,B,E are system matricies from a continuous time system,
% delta_t is the time step for discretization
% A1,B1,E1 are discretized matricies
% A = E1\A1
% B = E1\B1

% Math taken from "A collection of Benchmark examples for model reduction
% of linear time invariant dynamical systems" by Chahlaoui and Van Dooren

I = eye(length(A));

E1 = I-(delta_t/2)*A;
A1 = I+(delta_t/2)*A;
B1 = (delta_t/2)*B;

A = E1\A1;
B = (E1\B1);