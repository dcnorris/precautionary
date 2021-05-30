%% Implement BOIN designs using CCD machinery
:- use_module(library(clpz)).
:- use_module(library(lambda)).
:- use_module(library(precautionary/qbeta)).
:- use_module(library(precautionary/ccd)).

/*

1. Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials.
   J R Stat Soc C. 2015;64(3):507-523. doi:10.1111/rssc.12089

  Table 1 from [1]:

             cumulative patients treated at current dose (n_j):
              1    2    3    4    5    6    7    8    9    10    11    12

lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12

lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12

elimination   -    -   3/3  3/4  3/5  4/6  4/7  4/8  5/9  5/10  6/11  6/12

*/
% My initial emphasis is on generating all possible paths (CPE) for the BOIN
% design set forth in the table above. Although the BOIN design of [1] lacks
% any terminating principle except elimination of all doses, we do need such
% a rule here. The most natural type of rule, in view of the Table above,
% might be a 'stop-for-consensus' type of rule as found in package 'dtpcrm'.
% This is specified as a maximum number of patients to enroll at any 1 dose.

/*

A crucial task for this module is to generate the defining boundaries
for a whole range of BOIN designs. Without access to floating point math,
this involves essentially hard-coding each possible design. An analysis
of the design space thus proves essential.

This also should promote thought about how in general to represent such
design-defining boundaries, and communicate them to the CCD code.

Of course, representing trial designs as Prolog terms (who woulda thot?)
does present itself as one obvious solution. Ideally, this could be done
in such a way that ALL POSSIBLE CCDs could be fairly enumerated.
 
*/

/*

For convenience, we implement the simplest possible BOIN design described
in [1], namely the 'local BOIN' with priors 𝜋_0j = 𝜋_1j = 𝜋_2j, for which
Theorem 2 shows a long-term memory coherence property. Furthermore, we
adopt the recommended values 𝜙1 = 0.6𝜙, and 𝜙2 = 1.4𝜙.

*/

boin_targetpct_nmax(ccd(Elim, Deesc, Esc, Nmax), TargetPct, Nmax) :-
    phipct_lambda1_lambda2(TargetPct, Lambda1, Lambda2),
    findall(E, slope_floor(Lambda1, E), Esc_),
    findall(D, slope_ceiling(Lambda2, D), Deesc_),
    findall(R, elim_boundary(R, TargetPct), Elim_),
    ceiling_canonical(Elim_, Elim),
    ceiling_canonical(Deesc_, Deesc),
    floor_canonical(Esc_, Esc).

%?- boin_targetpct_nmax(BOIN, 25, 12).
%@    BOIN = ccd([3/5,4/8,5/10,6/12],[1/3,2/6,3/10,4/12],[0/1,1/6,2/11],12).

boin_targetpct_d_matrix(TargetPct, D, Matrix) :-
    Nmax = 12, % TODO: Somehow don't hard-code this
    boin_targetpct_nmax(BOIN, TargetPct, Nmax), % BOIN = ccd(...)
    ccd_d_matrix(BOIN, D, Matrix).

%?- Matrix+\(boin_targetpct_d_matrix(25, 3, Matrix)).
%@    Matrix = [[0/1,0/1]^[0/6]^[]~>3]
%@ ;  Matrix = [[0/1,0/1]^[1/6]^[]~>3]
%@ ;  Matrix = [[0/1,0/1]^[1/6]^[]~>3]
%@ ;  Matrix = [[0/1]^[0/2,2/6]^[]~>2]
%@ ;  Matrix = [[0/3]^[1/6,2/6]^[]~>2]
%@ ;  ...

%?- J+\(boin_targetpct_nmax(BOIN, 25, 12), D=1, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 0.092 seconds
%@    % CPU time: 0.096 seconds
%@    J = 10.

%?- J+\(boin_targetpct_nmax(BOIN, 25, 12), D=2, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 1.612 seconds
%@    % CPU time: 1.616 seconds
%@    J = 170.

%?- J+\(boin_targetpct_nmax(BOIN, 25, 12), D=3, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 9.009 seconds
%@    % CPU time: 9.014 seconds
%@    J = 949.

%?- J+\(boin_targetpct_nmax(BOIN, 25, 12), D=4, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 68.183 seconds
%@    % CPU time: 68.188 seconds
%@    J = 7139.

%?- J+\(boin_targetpct_nmax(BOIN, 25, 12), D=5, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 451.505 seconds
%@    % CPU time: 451.510 seconds
%@    J = 45927.

%% These are very close rational approximations to Eqs (2) & (3) in [1],
%% obtained using R function 'MASS:fractions' ...

phipct_lambda1_lambda2(15, 130904/1111271, 5280/29549).

phipct_lambda1_lambda2(20, 723/4598, 224107/939800).

phipct_lambda1_lambda2(25, 9043/45950, 33795/113257).

phipct_lambda1_lambda2(30, 18648/78853, 1172/3269).

%% From these 'slopes', I need to obtain floors & ceilings.
slope_floor(Y/X, T/N) :-
    N in 1..12, indomain(N),
    T #= (N * Y) // X.

slope_ceiling(Y/X, T/N) :-
    N in 1..12, indomain(N),
    T #= 1 + (N * Y) // X.

%?- slope_floor(1/3, T/6).
%@    T = 2.

%?- slope_ceiling(1/3, T/6).
%@    T = 3.

%?- findall(Q, slope_floor(1/3, Q), Floor), portray_clause(Floor).
%@ [0/1,0/2,1/3,1/4,1/5,2/6,2/7,2/8,3/9,3/10,3/11,4/12].
%@    Floor = [0/1,0/2,1/3,1/4,1/5,2/6,2/7,2/8,3/9,... / ...|...].

%% We also need to tabulate the 5% quantiles of the Beta distribution,
%% from which we will compute the elimination thresholds of [1,p515].

%% TODO: What are the implications for soundness of such a predicate,
%%       based like this on tabulation of a transcendental function?
elim_boundary(T/N, PhiPct) :-
    N in 3..12, % NB: Liu & Yuan (2015) eliminate only for N ≥ 3
    indomain(N),
    T in 1..N,
    post05_tally(Qupper, T/N),
    ratless(PhiPct/100, Qupper),
    Tminus1 #= T - 1,
    post05_tally(Qlower, Tminus1/N),
    ratless(Qlower, PhiPct/100).

%?- elim_boundary(T/N, 25).
%@ caught: error(existence_error(procedure,elim_boundary/2),elim_boundary/2)
%@    T = 3, N = 3
%@ ;  T = 3, N = 4
%@ ;  T = 3, N = 5
%@ ;  T = 4, N = 6
%@ ;  T = 4, N = 7
%@ ;  T = 4, N = 8
%@ ;  T = 5, N = 9
%@ ;  T = 5, N = 10
%@ ;  T = 6, N = 11
%@ ;  T = 6, N = 12
%@ ;  false.

%% Less-than relation on rationals
ratless(X1/Y1, X2/Y2) :-
    X1 * Y2 #< X2 * Y1. % 'cross-multiply'

%% Greater-than relation on rationals
ratmore(X1/Y1, X2/Y2) :- ratless(X2/Y2, X1/Y1).

%?- findall(T/N, elim_boundary(T/N, 25), Bdy), portray_clause(Bdy).
%@ [3/3,3/4,3/5,4/6,4/7,4/8,5/9,5/10,6/11,6/12].
%@    Bdy = [3/3,3/4,3/5,4/6,4/7,4/8,5/9,5/10,6/11,... / ...].

%% The 5% quantile of posterior probability of toxicity,
%% after observing toxicity tally T/N
post05_tally(P, T/N) :-
    Alpha #= T + 1,
    Beta #= N - T + 1,
    qbeta05_alpha_beta(P, Alpha, Beta).


%% From these predicates, we should now be able to obtain Prolog terms
%% that define the boundaries of a CCD.

%% TODO: Can I override library(ccd)'s (non-exported) definition of tally_decision/2?
%%  ANS: Apparently not! So it needs to be supplied as an argument.
/*
tally_decision(Q, Decision) :-
    tally_decision_ccd(Q, Decision, ccd([3/5, 4/8, 5/9, 6/12],
					[1/1, 2/4, 3/7, 4/8, 5/11],
					[0/1, 1/8],
					12)).
*/

%?- J+\(time(findall(Matrix, ccd_matrix(3, Matrix), _Paths)), length(_Paths, J)).
%@    % CPU time: 10.555 seconds
%@    % CPU time: 10.560 seconds
%@    J = 1151. % It seems my local tally_decision doesn't override; is this lexical scoping?
%@    % CPU time: 10.604 seconds
%@    % CPU time: 10.608 seconds
%@    J = 1151.
%@    % CPU time: 10.961 seconds
%@    % CPU time: 10.966 seconds
%@    J = 1151.
%@    % CPU time: 11.084 seconds
%@    % CPU time: 11.088 seconds
%@    J = 1151.
