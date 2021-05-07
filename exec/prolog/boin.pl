% [Adapted from previous 'rolling.pl', on the assumption
% that it represents my latest thinking on this subject,
% although it dealt specifically with extending 'aliquots'
% to the case of rolling enrollment.]

% Implementing BOIN within a larger class of designs.

/*

Having at last read [1] with an eye toward implementation, I see that BOIN
lends a formal status to a whole class of dose-escalation designs which may
prove extremely easy to implement and explore using Prolog.

Table 2 in the paper is especially helpful, suggesting that a general class
of transition rules may be relatively easy to specify in CLP terms, without
recourse to whatever Real-analytic computations might be used as heuristics
for *finding* them:

             cumulative patients treated at current dose (n_j):
              1    2    3    4    5    6    7    8    9    10    11    12

lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12

lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12

elimination   -    -   3/3  3/4  3/5  4/6  4/7  4/8  5/9  5/10  6/11  6/12


The dose-elimination concept also readily invites computations over *lists*,
which are subject to truncation in the same manner.

Furthermore, with the 'local BOIN' design having so few free parameters (indeed,
under the default choice delta = +/-40%, the target toxicity rate is all that we
have left!), we obtain the real possibility of caching precomputed T arrays for
a finite number of possible designs! Thus, we might cache BOIN T arrays on a grid
of (TTR,D) in {0.1, 0.15, 0.2, 0.25, 0.3} x {3, 4, 5, 6, 7, 8} which amounts to
only 30 distinct arrays. (An interesting question is whether it might suffice to
cache only the D=8 arrays, and then collapse these to the D<8 arrays by a quick
array computation. In this case, it might well be reasonable to cache on a finer
TTR grid.)

Also helpful for situating this type of design within a 'taxonomy' is [2], with
its notion of a 'cumulative cohort design'. This is a design where the escalation
decisions are made solely on the tally *at the current dose*, without regard to
the tallies at other doses. This is the basic pattern to be implemented here,
albeit with elimination of overly-toxic doses which departs from strict CCD form.


1. Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials.
   J R Stat Soc C. 2015;64(3):507-523. doi:10.1111/rssc.12089

2. Ivanova A, Flournoy N, Chung Y. Cumulative cohort design for dose-finding.
   Journal of Statistical Planning and Inference. 2007;137(7):2316-2327.
   doi:10.1016/j.jspi.2006.07.009

--------------------------------------------------------------------------------

*/
:- use_module(library(clpz)).
:- use_module(library(pio)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(lambda)).
:- use_module(library(time)).
%@    true.

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *).
*_.  % *Goal always succeeds

%% -----------------------------------------------------------------------------

% My initial emphasis is on generating all possible paths (CPE) for the BOIN
% design set forth in the table above. Although the BOIN design of [1] lacks
% any terminating principle except elimination of all doses, we do need such
% a rule here. The most natural formtype of rule, in view of the Table above,
% might be a 'stop-for-consensus' type of rule as found in package 'dtpcrm'.
% This is specified as a maximum number of patients to enroll at any 1 dose.

/*

Under a scheme of rolling enrollment, we no longer have the alternating sequence
of enrollment and assessment 'coroutines'. Rather, the GENERAL state is one where
EITHER further enrollment or resolution of pending DLT assessments may occur.

What if I plunged in, and explicitly used Prolog's distinction between ground and
uninstantiated variables? The sequence of enrolled patients is then a list of vars,
each in domain 0..1. The state transitions consist of either the instantiation of
non-ground vars, or else the extension of the list with new enrollment. One happy
consequence of this design is that it looks open to generalization to ordinal tox.

*/

/*

That looks promising! So the state representation now becomes the previous
left-right lists, except that the element of these lists are themselves lists!

*/

%% A 'dose cohort' (or 'cohort' for short) is a list of toxicity assessments
%% which are var upon enrollment, and become ground when completed.
%% (The association with a dose will be effected through the L^R stack pair
%% as previously implemented in the 'aliquots' design.)

cohort --> [].
cohort --> [T], { T in 0..1 }, cohort.

%% We will judge a cohort according to its *tally*, a proper fraction:
cohort_tally(C, T/N) :-
    length(C, N),
    C ins 0..1, % better than phrase(cohort, C)
    sum(C, #=, T). % clpz:sum/3 avoids awful "maplist(indomain, C), sum_list(C, T)"!

%?- phrase(cohort, C).
%@    C = []
%@ ;  C = [_A], clpz:(_A in 0..1)
%@ ;  C = [_A,_B], clpz:(_A in 0..1), clpz:(_B in 0..1)
%@ ;  C = [_A,_B,_C], clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1)
%@ ;  ...

%?- cohort_tally(C, T/N).
%@    C = [], T = 0, N = 0
%@ ;  C = [T], N = 1, clpz:(T in 0..1)
%@ ;  C = [_A,_B], N = 2, clpz:(_A+_B#=T), clpz:(T in 0..2), clpz:(_A in 0..1), clpz:(_B in 0..1)
%@ ;  C = [_A,_B,_C], N = 3, clpz:(_A+_B+_C#=T), clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1), clpz:(T in 0..3)
%@ ;  C = [_A,_B,_C,_D], N = 4, clpz:(_A+_B+_C+_D#=T), clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1), clpz:(_D in 0..1), clpz:(T in 0..4)
%@ ;  ...

%% Let's also be able to obtain best- and worst-case tallies
%% from a partially complete cohort of DLT assessments.
cohort_bestcase([], []).
cohort_bestcase([T|Ts], [B|Bs]) :-
    (	var(T),
	B #= 0
    ;	ground(T),
	B = T
    ),
    cohort_bestcase(Ts, Bs).

cohort_worstcase([], []).
cohort_worstcase([T|Ts], [W|Ws]) :-
    (	var(T),
	W #= 1
    ;	ground(T),
	W = T
    ),
    cohort_worstcase(Ts, Ws).

%?- cohort_bestcase([0,1,T], B).
%@    B = [0,1,0]
%@ ;  false.

%?- cohort_worstcase([0,1,T], W).
%@    W = [0,1,1]
%@ ;  false.

cohort_mintally(C, T/N) :-
    cohort_bestcase(C, B),
    cohort_tally(B, T/N).

cohort_maxtally(C, T/N) :-
    cohort_worstcase(C, W),
    cohort_tally(W, T/N).

%?- cohort_mintally([0,1,T], Q).
%@    Q = 1/3
%@ ;  false.

%?- cohort_maxtally([0,1,T], Q).
%@    Q = 2/3
%@ ;  false.

%% A fair enumeration of tallies will be helpful for developing
%% and testing comparison relations between tallies.

%% In general, I may require a max-enrollment *parameter*,
%% however unsightly it may be tagging along like this ...
tally_maxn(DLTs/Enrolled, MaxN) :-
    nonvar(MaxN),
    Enrolled in 0..MaxN,
    indomain(Enrolled),
    DLTs in 0..Enrolled,
    indomain(DLTs).

%?- tally_maxn(T/N, 6).
%@    T = 0, N = 0
%@ ;  T = 0, N = 1
%@ ;  T = 1, N = 1
%% ... as expected ...
%@ ;  T = 4, N = 6
%@ ;  T = 5, N = 6
%@ ;  T = 6, N = 6.

%% Having a 'default' max cohort size will make dev & test a bit easier:
tally(DLTs/Enrolled) :- tally_maxn(DLTs/Enrolled, 12).

%% Proof that tally/1 terminates!
%?- tally(_), false.
%@ false.

%?- tally(T/N).
%@    T = 0, N = 0
%@ ;  T = 0, N = 1
%@ ;  T = 1, N = 1
%@ ;  T = 0, N = 2
%@ ;  T = 1, N = 2
%% ... as expected ...
%@ ;  T = 8, N = 12
%@ ;  T = 9, N = 12
%@ ;  T = 10, N = 12
%@ ;  T = 11, N = 12
%@ ;  T = 12, N = 12.

%% How do tallies COMPARE?

/*
The essence of tally comparisons is one of REACHABILITY, i.e. EXISTENCE OF PATHS.
Generally, in fact, it is NON-REACHABILITY that supports any definite statement.
This might be stated correctly in terms of the DCG cohort//0. But this would not
be efficient. (Still, I could use such a statement as a correctness check!)

Strict inequalities like safer_than/2 mean that whichever argument is *earlier*
(i.e., has the smaller denominator) cannot under any conceivable path cross over
the other argument. Thus, the earlier is safer than the later if even a path of
100% toxicity that brings its denominator to the later's leaves its numerator
strictly less. Conversely, the later is safer if even a 100% non-toxic path from
the earlier does not make it look better than the later.
*/

safer_than(T0/N0, T1/N1) :-
    tally(T0/N0),
    tally(T1/N1),
    (	N0 #>= N1, T0 #< T1
    ;	N0 #< N1,
	MaxTox = N1 - N0,
	T0 + MaxTox #< T1
    ).

%?- safer_than(Q, 2/3).
%@    Q = 0/2
%@ ;  Q = 0/3
%@ ;  Q = 1/3
%@ ;  Q = 0/4
%@ ;  Q = 1/4
%@ ;  Q = 0/5
%@ ;  Q = 1/5
%@ ;  Q = 0/6
%@ ;  Q = 1/6
%@ ;  Q = 0/7
%@ ;  Q = 1/7
%@ ;  Q = 0/8
%@ ;  Q = 1/8
%@ ;  Q = 0/9
%@ ;  Q = 1/9
%@ ;  Q = 0/10
%@ ;  Q = 1/10
%@ ;  Q = 0/11
%@ ;  Q = 1/11
%@ ;  Q = 0/12
%@ ;  Q = 1/12
%@ ;  false.

%% Note e.g. that we cannot say safer_than(1/2, 2/3)
%% because we might enroll 1 with 1/2 ~~> 2/3, which
%% is not *strictly* less than 2/3.

:- op(900, xfx, &<).
&<(Q1, Q2) :- safer_than(Q1, Q2).
%?- Q &< 1/3.
%@    Q = 0/3
%@ ;  Q = 0/4
%@ ;  Q = 0/5
%@ ;  Q = 0/6
%@ ;  Q = 0/7
%@ ;  Q = 0/8
%@ ;  Q = 0/9
%@ ;  Q = 0/10
%@ ;  Q = 0/11
%@ ;  Q = 0/12
%@ ;  false.

%% Note that, under the new REACHABILITY semantics of tally comparisons,
%% we no longer make a claim such as 1/6 &< 1/3!
%% TODO: It seems likely that reachability semantics are more stringent
%%       than my previous formalism, so that we now can say fewer things.
%%       I ought to prove this---or let Prolog prove it for me.
%% TODO: What is the relationship (if any) between meta-interpretation
%%       and this change of tally-comparison semantics? Does MI affect
%%       only the IMPLEMENTATION---and not the meaning---of a program?
%% TODO: Given how far this reachability semantics departs from normal
%%       pharmacologic intuitions, perhaps even words like 'safer_than'
%%       should be reconsidered. I may well have abandoned all hope of
%%       capturing such intuitions in these predicates!

noworse_than(T0/N0, T1/N1) :-
    tally(T0/N0),
    tally(T1/N1),
    (	N0 #>= N1, T0 #=< T1
    ;	N0 #< N1,
	MaxTox = N1 - N0,
	T0 + MaxTox #=< T1
    ).

%?- safer_than(1/3, 2/3).
%@    true
%@ ;  false.

%?- noworse_than(1/3, 2/3).
%@    true
%@ ;  false.

%?- tally(1/3), tally(2/3), 3 #>= 3, 1 #=< 2.
%@    true.


%?- noworse_than(Q, 2/3).
%@    Q = 0/1
%@ ;  Q = 0/2
%@ ;  Q = 1/2
%@ ;  Q = 0/3
%@ ;  Q = 1/3
%@ ;  Q = 2/3
%@ ;  Q = 0/4
%@ ;  Q = 1/4
%@ ;  Q = 2/4
%@ ;  Q = 0/5
%@ ;  Q = 1/5
%@ ;  Q = 2/5
%@ ;  Q = 0/6
%@ ;  Q = 1/6
%@ ;  Q = 2/6
%@ ;  Q = 0/7
%@ ;  Q = 1/7
%@ ;  Q = 2/7
%@ ;  Q = 0/8
%@ ;  Q = 1/8
%@ ;  Q = 2/8
%@ ;  Q = 0/9
%@ ;  Q = 1/9
%@ ;  Q = 2/9
%@ ;  Q = 0/10
%@ ;  Q = 1/10
%@ ;  Q = 2/10
%@ ;  Q = 0/11
%@ ;  Q = 1/11
%@ ;  Q = 2/11
%@ ;  Q = 0/12
%@ ;  Q = 1/12
%@ ;  Q = 2/12
%@ ;  false.

%?- safer_than(T/N, 1/2).
%@    T = 0, N = 2
%@ ;  T = 0, N = 3
%@ ;  T = 0, N = 4
%@ ;  T = 0, N = 5
%@ ;  T = 0, N = 6
%@ ;  T = 0, N = 7
%@ ;  T = 0, N = 8
%@ ;  T = 0, N = 9
%@ ;  T = 0, N = 10
%@ ;  T = 0, N = 11
%@ ;  T = 0, N = 12
%@ ;  false.

:- op(900, xfx, &=<).
&=<(Q1, Q2) :- noworse_than(Q1, Q2).
%?- 1/3 &=< Q.
%@    Q = 1/1
%@ ;  Q = 1/2
%@ ;  Q = 2/2
%@ ;  Q = 1/3
%@ ;  Q = 2/3
%@ ;  Q = 3/3
%@ ;  Q = 2/4
%@ ;  Q = 3/4
%@ ;  Q = 4/4
%@ ;  Q = 3/5
%@ ;  Q = 4/5
%@ ;  Q = 5/5
%@ ;  Q = 4/6
%@ ;  Q = 5/6
%@ ;  Q = 6/6.

on_par(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0,
    tally(T1/N1), N1 #> 0,
    T0*N1 #= N0*T1.

:- op(900, xfx, &=).
&=(Q1, Q2) :- on_par(Q1, Q2).

%?- Q &= R.
%@    Q = 0/1, R = 0/1
%@ ;  Q = 0/1, R = 0/2
%@ ;  Q = 0/1, R = 0/3
%@ ;  Q = 0/1, R = 0/4
%@ ;  Q = 0/1, R = 0/5
%@ ;  Q = 0/1, R = 0/6
%@ ;  Q = 0/1, R = 0/7
%@ ;  Q = 0/1, R = 0/8
%@ ;  Q = 0/1, R = 0/9
%@ ;  Q = 0/1, R = 0/10
%@ ;  Q = 0/1, R = 0/11
%@ ;  Q = 0/1, R = 0/12
%@ ;  Q = 1/1, R = 1/1
%@ ;  Q = 1/1, R = 2/2
%@ ;  Q = 1/1, R = 3/3
%@ ;  Q = 1/1, R = 4/4
%@ ;  Q = 1/1, R = 5/5
%@ ;  Q = 1/1, R = 6/6
%@ ;  ...

%?- Q &< R.
%@    Q = 0/1, R = 1/1
%@ ;  Q = 0/1, R = 2/2
%@ ;  Q = 0/1, R = 3/3
%@ ;  Q = 0/1, R = 4/4
%@ ;  Q = 0/1, R = 5/5
%@ ;  Q = 0/1, R = 6/6
%@ ;  ...

% Strict tally-safer (&<) excludes tally-equals (&=)
%?- time((Q &< R, Q &= R)).
%@    % CPU time: 56.976 seconds
%@ false.

% (&=<) and '&>' are MUTUALLY EXCLUSIVE:
%?- time((Q &< R, R &=< Q)).
%@    % CPU time: 63.436 seconds
%@ false.


%% Here was my motivation for writing down on_par/2 in the first place:
%?- Q &= 1/3. %% 1/3 is (roughly) the 3+3 design's 'target toxicity rate'
%@    Q = 1/3
%@ ;  Q = 2/6
%@ ;  Q = 3/9
%@ ;  Q = 4/12
%@ ;  false.

%% TODO: Consider deriving all tally comparisons from cross-product zcompare,
%%       taking care to exclude 'incomparables' somehow.

%?- cohort_tally([0,0,T], Q), safer_than(Q, 2/3).
%@    Q = 0/3, T = 0
%@ ;  Q = 1/3, T = 1
%@ ;  false.

%?- Vs = [A,B,C], Vs ins 0..1, sum(Vs, #=, Sum), indomain(Sum).
%@    Vs = [0,A,A], A = 0, B = 0, Sum = 0, C = 0
%@ ;  Vs = [A,B,C], Sum = 1, clpz:(A+B+C#=1), clpz:(A in 0..1), clpz:(B in 0..1), clpz:(C in 0..1)
%@ ;  Vs = [A,B,C], Sum = 2, clpz:(A+B+C#=2), clpz:(A in 0..1), clpz:(B in 0..1), clpz:(C in 0..1)
%@ ;  Vs = [1,A,A], A = 1, B = 1, Sum = 3, C = 1.

%?- Vs = [A,B,C], Vs ins 0..1, sum(Vs, #=, Sum), labeling([], Vs).
%@    Vs = [0,A,A], A = 0, B = 0, Sum = 0, C = 0
%@ ;  Vs = [0,A,1], A = 0, B = 0, Sum = 1, C = 1
%@ ;  Vs = [0,1,A], A = 0, B = 1, Sum = 1, C = 0
%@ ;  Vs = [0,1,B], A = 0, B = 1, Sum = 2, C = 1
%@ ;  Vs = [1,0,B], A = 1, B = 0, Sum = 1, C = 0
%@ ;  Vs = [1,0,A], A = 1, B = 0, Sum = 2, C = 1
%@ ;  Vs = [1,A,0], A = 1, B = 1, Sum = 2, C = 0
%@ ;  Vs = [1,A,A], A = 1, B = 1, Sum = 3, C = 1.

%?- X = a, nonvar(X).
%@    X = a. 

%?- nonvar(X).
%@ false.

%% Therefore var/nonvar makes no sense declaratively!
%% var/1 and nonvar/1 are metalogical

%% So instead we can do this:
% Vs = [A,B,C], Vs ins 0..1, sum(Vs, #=, Sum), labeling([min(Sum)], Vs).
%% NB: If I want the first minimum, wrap this in once/1.

%% -----------------------------------------------------------------------------

% The BOIN rules suggest that comparisons of tallies might be generalized to
% comparisons between a tally and a *list* of tallies. Thus, e.g. the rule:
%
% lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12
%
% generates interest in comparisons such as T/N &=< [0/1, 0/2, 0/3, 0/4, ...].
%
% Does this comparison have an unambiguous or obvious meaning that might allow
% us to condense the list? It would be a shame to have to completely list all
% possible fractions in a sequence like this!
%
% With &=< at least, perhaps we would like to write instead:
%   T/N &=< [0/1, 1/8, 2/15],
% or possibly even abbreviate this as:
%   T/N &=< [1, 8, 15].
%
% The comparison would then hold if ANY of the mapped comparisons holds.
%
% What about the flip-side?
%
% lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12
%
% Here, we require comparisons like T/N &>= [1/1, 2/2, 2/3, ...] or possibly
% T/N &>= [1/1, 2/4, 3/6, ...] or even T/N &>= [1, 4, 6, ...].
%
% Again, this complex comparison holds if ANY component comparison holds.

% Define what a BOIN-type threshold is ...

zip([], [], []).
zip([X|Xs], [Y|Ys], [X/Y | XsYs]) :-
    zip(Xs, Ys, XsYs).

%?- zip([1,2,3], [4,5,6], Z).
%@    Z = [1/4,2/5,3/6].

%?- zip(T, N, [1/4,2/5,3/6]).
%@    T = [1,2,3], N = [4,5,6]
%@ ;  false.

tox_boundary(Bdy) :-
    zip(Ts, Ns, Bdy),
    sort(Ns, Ns), % TODO: Do I gain anything by thus removing degeneracy?
    maplist(\X^(X #> 0), Ns), % TODO: Is there a clpz idiom for this?
    %% Note there are further constraints on what might be considered
    %% reasonable tox boundaries. There shouldn't be duplicated Ns, e.g.,
    %% and all the Ts should be non-negative. But perhaps these are best
    %% postponed until some kind of *search* over the space of possible
    %% tox boundaries becomes desirable.
    maplist(#=<, Ts, Ns).

%?- tox_boundary([0/2, 1/4]).
%@    true
%@ ;  false.

%?- tox_boundary([0/5, 1/4]).
%@ false.

%% TODO: Why does sortedness of Ns not appear in constraints below?
%?- tox_boundary(Bdy).
%@    Bdy = []
%@ ;  Bdy = [_B/_A], clpz:(_A#>=_B), clpz:(_A in 1..sup)
%@ ;  Bdy = [_B/_A,_D/_C], clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_A in 1..sup), clpz:(_C in 1..sup)
%@ ;  Bdy = [_B/_A,_D/_C,_F/_E], clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_E#>=_F), clpz:(_A in 1..sup), clpz:(_C in 1..sup), clpz:(_E in 1..sup)
%@ ;  ...

/*
The question about a cumulative-cohort tally vis-à-vis a given threshold,
is whether it has touched or crossed the threshold toward 'extreme' values.
The past participle 'hit' seems reasonable, since it covers both the case
of touching and of having 'broken through'.
*/

%% TODO: Consider not leaning on (&<)-type tally comparisons, but rather
%%       going directly to CLP(Z) constraints. This will be faster, and
%%       also enable me to avoid superfluous choice points, etc.
%% OOH! Maybe I'm discovering that all the above efforts were unnecessary!
%%      Could it be that this file can just about start HERE?

/*
If I do 'cut out the middleman' in this way, what is the right representation
for the toxicity boundaries?

lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12
lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12
elimination   -    -   3/3  3/4  3/5  4/6  4/7  4/8  5/9  5/10  6/11  6/12

RELATIVE to a maximum cohort size of 12, the above boundaries gain sufficient
and minimal representation as the following lists:

lambda_1 [0/1, 1/8]
lambda_2 [1/1, 2/4, 3/6, 4/8, 5/10, 6/12]
elim     [3/5, 4/8, 5/10, 6/12]

This is because lambda_1 is a lower bounds, and so the removed limits were
non-binding. Likewise, since lambda_2 and elim were upper bounds the removed
limits were non-binding.

In general, minimal representations of tox boundaries will consist of lists
of tallies BdyRep satisfying:

minimal_maxn(BdyRep, MaxN) :-
    zip(Ts, Ns, BdyRep),
    T1 #=< Tn,
    Ts in T1..Tn,
    sort(Ns, Ns),
    maplist(\X^(X #=< MaxN), Ns).


*/

%% I unapologetically assume that the Bdy list is a minimal representation.
hit_floor(_/_, []) :- false.
hit_floor(T/N, [Tb/Nb | Bdys]) :-
    %minimal_maxn([Tb/Nb|Bdys], 12), %% TODO: Try activating this assertion
    %% TODO: The correctness of this predicate may depend only on the
    %%       sortedness of Ts and Ns. Work out the details!
    (	T #= Tb, N #>= Nb
    ;	T #> Tb, hit_floor(T/N, Bdys)
    ).
    
%?- hit_floor(0/3, [0/1, 0/2, 0/3, 0/4, 0/5, 0/6, 0/7, 1/8, 1/9, 1/10, 1/11, 1/12]).
%@    true
%@ ;  false.

%?- hit_floor(0/3, [0/1, 1/8]).
%@    true
%@ ;  false.

hit_ceiling(_/_, []) :- false.
hit_ceiling(T/N, [Tb/Nb | Bdys]) :-
    %minimal_maxn([Tb/Nb|Bdys], 12), %% TODO: Try activating this assertion
    (	T #= Tb, N #=< Nb
    ;	T #> Tb, hit_ceiling(T/N, Bdys)
    ).

%?- hit_ceiling(3/5, [1/1, 2/4, 3/6, 4/8, 5/10, 6/12]).
%@    true
%@ ;  false.


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Pick up from here...

%% Opportunity to replace the hard-coding here with meta-interpretation.
%% You might interpret 'more fairly', e.g.

%% ENROLL new patient 'T' into the lowest right-hand dose
state0_action_state(Ls ^ [R | Rs], enroll, Ls ^ [[T|R] | Rs]) :-
    enrollable_cohort(R),
    T in 0..1, % var T is pending tox assessment of the newly-enrolled patient
    %% TODO: impose restriction on acceptable tally
    true.

enrollable_cohort(R) :-
    length(R, N), N #< 6, % limit max cohort size
    true.


% Which tallies do not rule out further enrollment at a dose?
enrollable_tally(T0/N0) :-
    tally(T0/N0),
    (	N0 #= 0
    ;	N0 #= 3, T0 #=< 1
    ).
%?- enrollable_tally(Q).
%@ Q = 0/0 ;
%@ Q = _3234/3,
%@ _3234 in 0..1 ;
%@ false.

% Which tallies RULE OUT further enrollment at a dose?
unenrollable_tally(T/3) :- tally(T/3), T #> 1. % 'too toxic'
unenrollable_tally(T/6) :- tally(T/6). % 'enough already'

%?- state0_action_state(S0, enroll, S).
%@ S0 = _7610:[0/0|_7618],        % We enroll a new dose 0/0,
%@ S = _7610^[_7652/3|_7618],
%@ _7652 in 0..3 ;
%@ S0 = _11382:[_11394/3|_11390], % or an already-tried dose ..
%@ S = _11382^[_11424/6|_11390],
%@ _11394 in 0..1,                % of 0/3 or 1/3
%@ _11394+_11474#=_11424,         % and add correctly!
%@ _11474 in 0..3,                % This is # of new toxicities,
%@ _11424 in 0..4 ;               % which can sum to 4/6 at most.
%@ false.

/*
Instead of regarding merely TALLIES as presumably safe/toxic,
we apply these JUDGMENTS rather (and more properly) to the DOSES
themselves -- as represented through lists.
*/
presumably_safe([]). % [] corresponds to zero dose
presumably_safe([Q|_]) :- % "no need to enroll further at this dose to say it's okay"
    tally(Q),
    Q = T/6,   % TODO: Might I obtain the no-deescalation variant of 3+3
    T in 0..1. %       simply by modifying these 2 goals?

presumably_toxic([]) :- false.
presumably_toxic([Q|_]) :-
    tally(Q),
    Q = T/_,
    T #> 1.

%?- presumably_toxic([Q|_]).
%@ Q = _5352/3,
%@ _5352 in 2..3 ;
%@ Q = _6862/6,
%@ _6862 in 2..6.

%% PROOF: Presumable tallies are NOT enrollable
%?- (presumably_safe([T]); presumably_toxic([T])), enrollable_tally(T).
%@ false.

%% PROOF: Enrollable and unenrollable tallies are disjoint
%?- enrollable_tally(Q), unenrollable_tally(Q).
%@ false.

%?- enrollable_tally(T/N).
%@ T = N, N = 0 ;
%@ N = 3,
%@ T in 0..1 ;
%@ false.

% Handling the (exceptional) MTD-not-found case in the ENROLLMENT
% phase is formally quite natural, even if it appears odd from a
% perspective that ':' is the state upon which the 'routine' trial
% operations act. (The literal interpretation of ':' as the phase
% where trial administrative staff are at work is thus misleading.)
%% AHA! This demands modification along with the 'conspicuous'
%% clause below.
%% state0_action_state(Ls : [], stop, mtd_notfound(D)) :-
%%     presumably_safe(Ls),
%%     length(Ls, D). % RP2D is highest prespecified dose, D.
% This is by far the MOST CONSPICUOUS clause of state0_action_state/3,
% since it clearly catches an edge case _:[] and changes the stacks.
% There is a real sense in which this exceptional case is handled as
% an 'oopsie' here, rather than straightforwardly.
% Perhaps we see here a tension between 'elegance' and 'bluntness',
% BOTH of which seem like attrbutes of good Prolog programming.
% ---
% Now that I find this clause causing problems during DCG translation
% to the compact D{^,-,*}T notation, I believe a special term is
% warranted to distinguish this case.
%% state0_action_state([Q0|Ls] : [], enroll, Ls ^ [Q]) :-
%%     tally0_cohort_tally(Q0, _, Q).


%% JUDGE the toxicities, and WEIGH the dose-escalation decision.
%% NB: We require a dose R that we can escalate TO!
state0_action_state(Ls ^ [Q,R|Rs], escalate, [Q|Ls] : [R|Rs]) :-
    Q &< 1/3.
% In the case where we WOULD ESCALATE, BUT CAN'T, we record this
% as a 'clamp' action, by analogy with voltage clamping.
state0_action_state(Ls ^ [Q], clamp, Ls : [Q]) :-
    Q &< 1/3,
    enrollable_tally(Q).
state0_action_state(Ls^[Q|Rs], stay, Ls:[Q|Rs]) :-
    Q &= 1/3, % <-- Roughly a statement about meeting 'target toxicity rate'
    enrollable_tally(Q).
state0_action_state([L|Ls] ^ Rs, deescalate, Ls : [L]) :-
    enrollable_tally(L),
    presumably_toxic(Rs).
state0_action_state(Ls ^ [R|Rs], stop, declare_mtd(MTD)) :-
    presumably_safe(Ls),
    presumably_toxic([R|Rs]),
    length(Ls, MTD).
% Here again we treat specially the case where 1st arg is Ls ^ [Q] ...
state0_action_state(Ls ^ [Q], stop, declare_mtd(MTD)) :-
    presumably_safe([Q|Ls]),
    length([Q|Ls], MTD).


% What does a  cohort look like?
cohort(DLTs/N) :-
    N #= 3, % initially, we permit only cohorts of 3
    DLTs in 0..N,
    *indomain(DLTs).
%?- cohort(C).
%@ C = _44170/3,
%@ _44170 in 0..3.


%?- unenrollable_tally(Q).
%@ Q = _5464/3,
%@ _5464 in 2..sup ;
%@ Q = _6286/6.

% How do cohorts ACCUMULATE into tallies?
tally0_cohort_tally(T0/N0, T_/N_, T/N) :-
    enrollable_tally(T0/N0),
    cohort(T_/N_),
    tally(T/N),
    T #= T0 + T_,
    N #= N0 + N_.
%?- tally0_cohort_tally(T0, C, T).
%@ T0 = 0/0,
%@ C = T, T = _10336/3,
%@ _10336 in 0..3 ;
%@ T0 = _14210/3,
%@ C = _14228/3,
%@ T = _14246/6,
%@ _14210 in 0..1,
%@ _14210+_14228#=_14246,
%@ _14228 in 0..3,
%@ _14246 in 0..4 ;
%@ false.


/* A truly ESSENTIAL core pharmacologic concept undergirding dose-finding
 * trials is the monotonicity of the dose-toxicity function. To manifest
 * this notion in our code, we model the DOSE as a DESCENDING LIST.
 * (Cf. the treatment of the hyperreals as convergent sequences of reals.)
 * In this treatment, it is not necessary to carry along dose labels in
 * these lists. The mere STRUCTURE of the list suffices to enable integer
 * doses to be 'read off' via length/2. Accordingly, we are free to use
 * the list elements to maintain tallies of toxicities observed at the
 * several doses of the trial.
*/

%% NB: You can do this instead, since tally/1 terminates:
%% tallies(Qs) :- maplist(tally, Qs).

% Describe a LIST of TALLIES for consecutive prespecified doses:
tallies([]).
tallies([Q|Qs]) :-
    length(Qs, _), % enumerate solutions fairly
    tally(Q),
    tallies(Qs).
%?- length(C, 1), tallies(C).
%@ C = [0/0] ;
%@ C = [_7170/3],
%@ _7170 in 0..3 ;
%@ C = [_8490/6],
%@ _8490 in 0..6. % NB: Only 0..3 could occur in a 3+3 trial

%?- length(C, 2), tallies(C), false. % Does it terminate?
%@ false.                            % Yes.

/* An assumption we make about the dose-escalation trials modeled here
 * is that AT ANY GIVEN MOMENT they always partition the available doses
 * into 2 groups:
 * 'L' - A descending list of doses considered 'too Low' to enroll RIGHT NOW
 * 'R' - An ascending list of doses of which the head looks Right to enroll.
 * We can imagine these lists as 2 stacks, on the 'left' and 'right'.
 * 
 * Since it proves useful to describe the trial in two alternating phases,
 * we designate the state of the trial via L and R lists conjoined with
 * alternating infix functors:
 * '^' - Imagine: 'scales of Justice' where we make judgements about tallies
 * ':' - Imagine: an ellipsis (or dripping faucet) where we enroll patients
 *       (who 'trickle in'), wait for & assess toxicities.
 */



/****
The following general queries seem to PROVE
key properties of 3+3 have been attained.
 ****/
%?- RP2D in 0..sup, state0_action_state(S0, stop, mtd_notfound(RP2D)).
%@ RP2D = 0,
%@ S0 = []:[] ;
%@ RP2D = 1,
%@ S0 = [_6146/6]:[],
%@ _6146 in 0..1 ;
%@ RP2D = 2,
%@ S0 = [_7652/6, _7658]:[],
%@ _7652 in 0..1 ;
%@ RP2D = 3,
%@ S0 = [_9170/6, _9176, _9182]:[],
%@ _9170 in 0..1 ; ...

%?- MTD in 0..sup, state0_action_state(S0, stop, declare_mtd(MTD)).
%@ MTD = 0,
%@ S0 = []^[_2400/3|_2396],
%@ _2400 in 2..3 ;
%@ MTD = 0,
%@ S0 = []^[_4346/6|_4342],
%@ _4346 in 2..6 ;
%@ MTD = 1,
%@ S0 = [_7710/6]^[_7722/3|_7718],
%@ _7710 in 0..1,
%@ _7722 in 2..3 ;
%@ MTD = 2,
%@ S0 = [_9476/6, _9482]^[_9494/3|_9490],
%@ _9476 in 0..1,
%@ _9494 in 2..3 ;
%@ MTD = 3,
%@ S0 = [_11254/6, _11260, _11266]^[_11278/3|_11274],
%@ _11254 in 0..1,
%@ _11278 in 2..3 ; ...

%?- MTD in 0..3, indomain(MTD), state0_action_state(S0, stop, declare_mtd(MTD)).
%@ MTD = 0,
%@ S0 = []^[_2842/3|_2838],
%@ _2842 in 2..3 ;
%@ MTD = 0,
%@ S0 = []^[_4762/6|_4758],
%@ _4762 in 2..6 ;
%@ MTD = 1,
%@ S0 = [_8222/6]^[_8234/3|_8230],
%@ _8222 in 0..1,
%@ _8234 in 2..3 ;
%@ MTD = 1,
%@ S0 = [_10396/6]^[_10408/6|_10404],
%@ _10396 in 0..1,
%@ _10408 in 2..6 ;
%@ MTD = 2,
%@ S0 = [_14060/6, _14066]^[_14078/3|_14074],
%@ _14060 in 0..1,
%@ _14078 in 2..3 ;
%@ MTD = 2,
%@ S0 = [_16252/6, _16258]^[_16270/6|_16266],
%@ _16252 in 0..1,
%@ _16270 in 2..6 ;
%@ MTD = 3,
%@ S0 = [_19894/6, _19900, _19906]^[_19918/3|_19914],
%@ _19894 in 0..1,
%@ _19918 in 2..3 ;
%@ MTD = 3,
%@ S0 = [_22104/6, _22110, _22116]^[_22128/6|_22124],
%@ _22104 in 0..1,
%@ _22128 in 2..6.
/* ---
 This makes a very nice characterization of the 'typical' or 'motivating'
 case of successful conclusion of the trial. In this variant of 3+3, the
 most desirable conclusion is to have observed {0,1}/6 toxities at the
 declared MTD, and excessive toxicity (2 or more) at the next-higher dose.
 */

%?- state0_action_state(S0, stop, mtd_notfound(MTD)).
%@ S0 = []:[],
%@ MTD = 0 ;
%@ S0 = [_1510]:[],
%@ MTD = 1 ;
%@ S0 = [_1510, _2558]:[],
%@ MTD = 2 ;
%@ S0 = [_1510, _2558, _3606]:[],
%@ MTD = 3 .

%% DISCUSS: Non-terminating. From a verification perspective,
%% this is not so troubling, since it merely asks Prolog to
%% confirm something we can already verify BY INSPECTION.
%?- state0_action_state(_:[], A, _), dif(A, stop).
%@ Action (h for help) ? abort
%@ % Execution Aborted

%?- state0_action_state(S0, deescalate, S).
%@ S0 = [0/3|_9574]^[[_9596/3|_9592]|_9586],
%@ S = _9574:[0/3],
%@ _9596 in 2..3 ;
%@ S0 = [0/3|_11542]^[[_11564/6|_11560]|_11554],
%@ S = _11542:[0/3],
%@ _11564 in 2..6.
/* ---
 What does this show about de-escalation?
 1. It only ever occurs back to a dose where tally is 0/3
 2. It never occurs after any of {0,1}/3, {0,1}/6
 3. Resulting state S is of 'pending' type ":"
 4. This 0/3 tally is the 'focus dose' at head of right-hand list
 5. The dose we de-escalated FROM looked like T/6 with T in 3..6
 6. The dose we de-escalated FROM gets TRUNCATED from right-hand list
 */

%?- state0_action_state(S0, escalate, S).
%@ S0 = _504^[0/3|_512],
%@ S = [0/3|_504]:_512 ;
%@ S0 = _4778^[_4790/6|_4786],
%@ S = [_4790/6|_4778]:_4786,
%@ _4790 in 0..1 ;
%@ false.
/* ---
 What does this say about escalation?
 1. It may occur ONLY after a 0/3, 0/6 or 1/6 observation
    (NB: The 0/6 observation is not excluded BY THIS PREDICATE ALONE.)
 2. Escalation leaves us in a ':' state
 3. The RHS of the ':' state MAY be empty list [].
 */

%?- state0_action_state(S0, stay, S).
%@ S0 = _7196^[1/3|_7204],
%@ S = _7196:[1/3|_7204] ;
%@ false.
/* ---
 That result speaks for itself!
 */


%?- X is 4 rdiv 7.
%@ X = 4r7.

%?- X is log(3).
%@ X = 1.0986122886681098.
%@ ERROR: Syntax error: Operator expected
%@ ERROR: X is log
%@ ERROR: ** here **
%@ ERROR:  3 . 

%% About which (indeterminate) tallies do we presume nothing?
%?- tally(I), \+ presumably_safe([I]), \+ presumably_toxic([I]).
%@ I = 0/0 ;
%@ I = 0/3 ;
%@ I = 1/3 ;
%@ false.
%% NOTE: These are precisely the enrollable doses!

%?- enrollable_tally(T).
%@ T = 0/0 ;
%@ T = _8562/3,
%@ _8562 in 0..1 ;
%@ false.

