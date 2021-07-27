%% Implement BOIN designs using CCD machinery
:- use_module(library(clpz)).
:- use_module(library(lambda)).
:- use_module(library(format)).
:- use_module(qbeta).
:- use_module(benchmarking).

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
% rules here. Two rule are provided, implemented via the 4th & 5th arguments
% of the ccd(_, _, _, Cmax, Nmax) term.

:- use_module(ccd).

/*

For convenience, we implement the simplest possible BOIN design described
in [1], namely the 'local BOIN' with priors 𝜋_0j = 𝜋_1j = 𝜋_2j, for which
Theorem 2 shows a long-term memory coherence property. Furthermore, we
adopt the recommended values 𝜙1 = 0.6𝜙, and 𝜙2 = 1.4𝜙.

*/

boin_targetpct_cmax_nmax(ccd(Elim, Deesc, Esc, FullCoh, Nmax), TargetPct, FullCoh, Nmax) :-
    phipct_lambda1_lambda2(TargetPct, Lambda1, Lambda2),
    findall(E, slope_floor(Lambda1, E), Esc_),
    findall(D, slope_ceiling(Lambda2, D), Deesc_),
    findall(R, elim_boundary(R, TargetPct), Elim_),
    ceiling_canonical(Elim_, Elim),
    ceiling_canonical(Deesc_, Deesc),
    floor_canonical(Esc_, Esc).

%?- boin_targetpct_cmax_nmax(BOIN, 25, 6, 12).
%@    BOIN = ccd([3/5,4/8,5/10,6/12],[1/3,2/6,3/10,4/12],[0/1,1/6,2/11],6,12).

boin_targetpct_d_cmax_nmax_matrix(TargetPct, D, CohortMax, EnrollMax, Matrix) :-
    boin_targetpct_cmax_nmax(BOIN, TargetPct, CohortMax, EnrollMax),
    ccd_d_matrix(BOIN, D, Matrix).

%?- Matrix+\(boin_targetpct_d_cmax_nmax_matrix(25, 3, 6, 24, Matrix)).
%@    Matrix = ([0/1,0/1]^[0/6]^[]~>3)
%@ ;  Matrix = ([0/1,0/1]^[1/6]^[]~>3)
%@ ;  Matrix = ([0/1,0/1]^[1/6]^[]~>3)
%@ ;  Matrix = ([0/1]^[0/2,2/6]^[]~>2)
%@ ;  Matrix = ([0/3]^[1/6,2/6]^[]~>2)
%@ ;  ...

%?- J+\(boin_targetpct_cmax_nmax(BOIN, 25, 6, 12), D=1, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 0.619 seconds
%@    % CPU time: 0.623 seconds
%@    J = 10.

%?- J+\(boin_targetpct_cmax_nmax(BOIN, 25, 6, 12), D=2, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 11.246 seconds
%@    % CPU time: 11.250 seconds
%@    J = 170.
%@    % CPU time: 9.888 seconds
%@    % CPU time: 9.893 seconds
%@    J = 170.
%@    % CPU time: 9.745 seconds
%@    % CPU time: 9.749 seconds
%@    J = 170.
%@    % CPU time: 7.470 seconds
%@    % CPU time: 7.476 seconds
%@    J = 170
%@ ;  false.
%@    % CPU time: 5.684 seconds
%@    % CPU time: 5.689 seconds
%@    J = 170
%@ ;  false.
%@    % CPU time: 1.612 seconds
%@    % CPU time: 1.616 seconds
%@    J = 170.

%?- J+\(boin_targetpct_cmax_nmax(BOIN, 25, 6, 24), D=3, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 73.365 seconds
%@    % CPU time: 73.370 seconds
%@    J = 949.
%@    % CPU time: 60.834 seconds
%@    % CPU time: 60.838 seconds
%@    J = 949.

%?- J+\(boin_targetpct_cmax_nmax(BOIN, 25, 6, 24), D=4, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 600.243 seconds
%@    % CPU time: 600.247 seconds
%@    J = 7139.
%@    % CPU time: 396.386 seconds
%@    % CPU time: 396.390 seconds
%@    J = 7139
%@ ;  false.
%@    % CPU time: 231.354 seconds
%@    % CPU time: 231.359 seconds
%@    J = 7139
%@ ;  false.
%@    % CPU time: 68.183 seconds
%@    % CPU time: 68.188 seconds
%@    J = 7139.

%?- J+\(boin_targetpct_cmax_nmax(BOIN, 25, 6, 24), D=5, time(findall(M, ccd_d_matrix(BOIN, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 2574.457 seconds
%@    % CPU time: 2574.462 seconds
%@    J = 31475.
%@    % CPU time: 1822.763 seconds
%@    % CPU time: 1822.767 seconds
%@    J = 31475 % < 45927, as expected since now D*6 > 24
%@ ;  false.
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
    #T #= (#N * #Y) // #X.

slope_ceiling(Y/X, T/N) :-
    N in 1..12, indomain(N),
    #T #= 1 + (#N * #Y) // #X.

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
    #Tminus1 #= #T - 1,
    post05_tally(Qlower, Tminus1/N),
    ratless(Qlower, PhiPct/100).

%?- elim_boundary(T/N, 25).
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
    #X1 * #Y2 #< #X2 * #Y1. % 'cross-multiply'

%% Greater-than relation on rationals
ratmore(X1/Y1, X2/Y2) :- ratless(X2/Y2, X1/Y1).

%% The 5% quantile of posterior probability of toxicity,
%% after observing toxicity tally T/N
post05_tally(P, T/N) :-
    #Alpha #= #T + 1,
    #Beta #= #N - #T + 1,
    qbeta05_alpha_beta(P, Alpha, Beta).


%% Write out tab-delimited BOIN path matrices BOIN<TgtPct>_<D>_<CohortMax>_<EnrollMax>.tab

write_T(TargetPct, D, CohortMax, EnrollMax) :-
    phrase(format_("BOIN~d-~d-~d-~d.tab",
		   [TargetPct, D, CohortMax, EnrollMax]),
	   Filename),
    boin_targetpct_cmax_nmax(BOIN, TargetPct, CohortMax, EnrollMax),
    ccd_d_tabfile(BOIN, D, Filename).

%?- assertz(clpz:monotonic).
%@    true.

%?- write_T(25, 2, 6, 24).
%@ Opening file 'BOIN25-2-6-24.tab'...
%@ Writing path vectors .... done (0.31 minutes).
%@    true.
%@ Opening file 'BOIN25-2-6-24.tab'...
%@ Writing path vectors .... done (0.20 minutes). J = 170
%@    true.

%?- write_T(25, 3, 6, 24).
%@ Opening file 'BOIN25-3-6-24.tab'...
%@ Writing path vectors .... done (1.89 minutes).
%@    true.
%@ Opening file 'BOIN25-3-6-24.tab'...
%@ Writing path vectors .... done (1.22 minutes). J = 949
%@    true.

%?- write_T(25, 4, 6, 24).
%@ Opening file 'BOIN25-4-6-24.tab'...
%@ Writing path vectors .... done (16.29 minutes).
%@    true.
%@ Opening file 'BOIN25-4-6-24.tab'...
%@ Writing path vectors .... done (10.37 minutes).
%@    true.

%?- write_T(25, 4, 9, 24).
%@ Opening file 'BOIN25-4-9-24.tab'...
%@ Writing path vectors .... done (177.76 minutes). J = 131,898
%@    true.

% NB: Net rate is 700+ paths/minute.
