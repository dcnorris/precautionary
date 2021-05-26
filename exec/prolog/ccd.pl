% Abstracting the cumulative-cohort design (CCD) principle for dose escalation

/*

The 'cumulative-cohort design' (CCD) concept in [1] establishes a class of
dose-escalation designs in which escalation decisions are a function strictly
of the cumulative toxicity rate AT THE CURRENT DOSE. These designs are thus
positioned midway between what [1] calls 'group designs' (in which escalation
depends only on the rate observed in the current COHORT), and general designs
where these decisions depend on the observed rates AT ALL DOSES. (Note that,
in the dose-escalation context where patients are treated as *exchangeable*,
the list of dose-wise tallies [T1/N1, T2/N2, ..., Td/Nd] indeed a sufficient
statistic for the full history of the trial.)

With its 'middling' size, this CCD class seems likely to be amenable to modes
of analysis that Prolog and CLP(ℤ) can bring to bear. But also apparently CCD
is general enough to encompass BOIN [2] and other interval-based methods such
as mTPI [3].

Indeed, even the analytical approach pursued in [1] invites comparison with
what we may attempt here by quite different---and, I think, more powerful &
suitable---methods. Whereas Ivanova, Flournoy & Chung attempt to reduce CCDs
to locally equivalent group designs, we might aim toward a similar analysis
of CRM and other such designs via their local approximation as CCDs!

Moreover, the effort in [1] towards gaining sharp control over the cardinality
of the design space also anticipates some of the challenges we will face here.

What we here call 'CCD' owes something to the generalizations drawn from [2],
especially the argument (starting at the bottom of p.514) that dose removal
boundaries are generally necessary in these designs.

1. Ivanova A, Flournoy N, Chung Y. Cumulative cohort design for dose-finding.
   Journal of Statistical Planning and Inference. 2007;137(7):2316-2327.
   doi:10.1016/j.jspi.2006.07.009

2. Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials.
   J R Stat Soc C. 2015;64(3):507-523. doi:10.1111/rssc.12089

3. Ji Y, Liu P, Li Y, Bekele BN. A modified toxicity probability interval method
   for dose-finding trials. Clinical Trials. 2010;7(6):653-663. doi:10.1177/1740774510382799

--------------------------------------------------------------------------------

While it seems appropriate to retain the apt term "cumulative cohort design",
we can go beyond the limited terms in which CCDs were defined and optimized
in [1]. The implementation of [2] actually points the way here.

Table 2 in Liu & Yuan is especially helpful, suggesting that a general class
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


--------------------------------------------------------------------------------

*/
:- use_module(library(clpz)).
:- use_module(library(pio)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(time)).
:- use_module(library(debug)).
:- use_module(library(dif)).
:- use_module(library(reif)).
:- use_module(library(format)).
:- use_module(library(lambda)).

% TODO: Review library(debug); it includes * and other useful predicates
% TODO: Look at clpz:automaton/_ predicates
% TODO: Look at clpz:parse_clpz, with its one-line-per-thing-that-can-occur.
%       Is there a similar language for cumulative-cohort designs?
%       Try invoking clpz:make_parse... clauses.

%?- clpz:make_parse_clpz(Clauses), maplist(portray_clause, Clauses).
%@ parse_clpz(A,B) :-
%@    cyclic_term(A),
%@    !,
%@    domain_error(clpz_expression,A).
%@ parse_clpz(A,B) :-
%@    var(A),
%% ...

%?- clpz:make_parse_reified(Clauses), maplist(portray_clause, Clauses).
%@ parse_reified_clpz(A,B,C,D,E) :-
%@    cyclic_term(A),
%@    D=F,
%@    !,
%@    F=G,
%@    domain_error(clpz_expression,A),
%@    G=H,
%@    a(C,H,E).
%@ parse_reified_clpz(A,B,C,D,E) :-
%@    var(A),
%% ...

% TODO: Look at the term_expansion stuff at the bottom of clpz.pl

%% -----------------------------------------------------------------------------

%% The most fundamental relation between tallies is REACHABILITY,
%% the question of whether a possible path exists connecting them.
qcompare(~~, T1/N1, T2/N2) :- % REACHABILITY
    %% The following conditions translate as follows, under 3 scenarios:
    %% (=) If N1==N2, then we get T2 =< T1 =< T2 which holds iff T1==T2.
    %% (<) If N1 < N2, these conds translate to T1 =< T2 =< T1 + MaxTox.
    %% (>) If N1 > N2, these conds translate to T2 =< T1 =< T2 + MaxTox.
    T2 #=< T1 + max(0, N2 - N1),
    T1 #=< T2 + max(0, N1 - N2).

/*

This ~~ relation naturally gives us EQUIVALENCE CLASSES that we
can use to define a PARTIAL ORDERING on tallies that extends the
natural ordering between 2 tallies sharing a common denominator.
Specifically, for any 2 tallies, Q1 = T1/N1 and Q2 = T2/N2 with
N1 < N2 (we call Q1 the 'earlier' tally and Q2 the 'later' one),
we can project the earlier tally forward in time to a reachable
of tallies {Tr/N2 | Tr/N2 ~~ T1/N1}, and then compare that set
(in the obvious way) to Q2. If every element of the set bears
some relation to Q2, then we can say Q1 bears that relation
modulo reachability.

A graphical analysis of this relation is helpful. Here are the
relations generated by the tally 3/5, on the simplex of tallies
depicted as a grid:

  7  >  >  >  >  >  >  >  >  >  ≥

  6  >  >  >  >  >  >  >  >  ≥

  5  >  >  >  >  >  >  >  ≥

  4  >  >  >  >  >  >  ≥

  3  ≥  ≥  ≥  ≥  ≥  =  ≤  ≤  ≤  ≤

  2              ≤  <  <  <  <  <

  1           ≤  <  <  <  <  <  < 

  0        ≤  <  <  <  <  <  <  <
 T
     0  1  2  3  4  5  6  7  8  9
    N

Note that nonexistent tallies (T/N with T>N) are included in the
relation, in order to emphasize the geometry, especially in the
'northwest' part of the figure, where imposing T≯N would wipe out
the ramp function formed by the (≥). Note how an (=)
is of course generated at 3/5, and how paired rays of (≤) and (≥)
emanate from 3/5, defining boundaries with interiors where the
corresponding *strict* relations hold. Note also that there are
two acute angles where neither relation holds. Thus, for example,
we cannot assert any relation between 2/3 and 3/5, because the
tallies reachable from 2/3 span the range 2/5--4/5, bracketed by
2/5 < 3/5 < 4/5. (Indeed, the reachability relation does hold,
so we could write 2/3 ~~ 3/5. But, as will become clear below,
it is the inequalities (≤) and (≥) that will provide the motive
force driving dose-escalation decisions. So reachability ~~ is
not emphasized in this figure.)

The geometry of this figure---and especially the CONVEXITY of the
obtuse-angle regions bounded by (≤) and (≥)---aids in understanding
the implementation of the qcompare/3 clauses below, and appreciating
properties of these relations, such as their transitivity.

*/

%% Note that, for reasons of performance, we express these relations
%% in terms of CLP(ℤ) constraints rather than by directly generating
%% the reachable sets.
%% TODO: Nevertheless, Prolog-based PROOFS via such representations
%%       might well have intrinsic interest.
%% Note also that we implement these comparisons on all of ℕ × ℕ, so that
%% queries about tallies must assert the simplex constraint themselves.
qcompare(=<, T1/N1, T2/N2) :-
    T1 + max(0, N2 - N1) #=< T2.
	
qcompare(<, T1/N1, T2/N2) :-
    T1 + max(0, N2 - N1) #< T2.
	
qcompare(>=, T1/N1, T2/N2) :-
    T1 #>= T2 + max(0, N1 - N2).

qcompare(>, T1/N1, T2/N2) :-
    T1 #> T2 + max(0, N1 - N2).

%% Reified versions of the above, as done at bottom of clpz.pl
qcompare(=<, T1/N1, T2/N2, Truth) :-
    T1 + max(0, N2 - N1) #=< T2 #<==> B,
    zo_t(B, Truth).
	
qcompare(<, T1/N1, T2/N2, Truth) :-
    T1 + max(0, N2 - N1) #< T2 #<==> B,
    zo_t(B, Truth).
	
qcompare(>=, T1/N1, T2/N2, Truth) :-
    T1 #>= T2 + max(0, N1 - N2) #<==> B,
    zo_t(B, Truth).

qcompare(>, T1/N1, T2/N2, Truth) :-
    T1 #> T2 + max(0, N1 - N2) #<==> B,
    zo_t(B, Truth).

zo_t(0, false).
zo_t(1, true).

%% Operators for the truly useful comparisons of BOIN:
:- op(900, xfx, &=<).
&=<(Q1, Q2) :-
    qcompare(=<, Q1, Q2).

&=<(Q1, Q2, Truth) :- % reified
    qcompare(=<, Q1, Q2, Truth).

:- op(900, xfx, &>=).
&>=(Q1, Q2) :-
    qcompare(>=, Q1, Q2).

&>=(Q1, Q2, Truth) :- % reified
    qcompare(>=, Q1, Q2, Truth).


%% Fairly enumerate tally pairs (Q1, Q2) pairs JOINTLY to avoid
%% repeating tests while 'sweeping' the space of cases to exclude
%% the 'inconceivable' systematically on increasing subsets of
%% the domain.
%% This is intended to be invoked with (Size in Min..Max), as e.g.
%% by inconceivable/3.
tally_pair(T1/N1, T2/N2, Size) :-
    Size #> 0,
    N1 #= Size, % NB: naive (N1 in 1..Size) would duplicate pairs
    N2 in 1..N1, indomain(N2),
    T1 in 0..N1,
    T2 in 0..N2.
    %%*indomain(T1), indomain(T2). % TODO: Find slicker way to toggle labeling

%?- Size in 1..3, indomain(Size), tally_pair(Q1, Q2, Size).
%@    Size = 1, Q1 = _A/1, Q2 = 0/1, clpz:(_A in 0..1)
%@ ;  Size = 1, Q1 = _A/1, Q2 = 1/Size, clpz:(_A in 0..1)
%@ ;  Size = 2, Q1 = _A/2, Q2 = 0/1, clpz:(_A in 0..2)
%@ ;  Size = 2, Q1 = _A/2, Q2 = 1/1, clpz:(_A in 0..2)
%@ ;  Size = 2, Q1 = _A/2, Q2 = 0/2, clpz:(_A in 0..2)
%@ ;  Size = 2, Q1 = _A/2, Q2 = 1/2, clpz:(_A in 0..2)
%@ ;  Size = 2, Q1 = _A/2, Q2 = 2/Size, clpz:(_A in 0..2)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 0/1, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 1/1, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 0/2, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 1/2, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 2/2, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 0/3, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 1/3, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 2/3, clpz:(_A in 0..3)
%@ ;  Size = 3, Q1 = _A/3, Q2 = 3/Size, clpz:(_A in 0..3).

%% Demonstrate that strict inequalities are exclusive of ~~
%?- inconceivable((tally_pair(Q1, Q2, MaxN), qcompare(>, Q1, Q2), qcompare(~~, Q1, Q2)), MaxN, 1..12).
%@  % MaxN = 1 ...   % CPU time: 0.020 seconds
%@  % MaxN = 2 ...   % CPU time: 0.062 seconds
%@  % MaxN = 3 ...   % CPU time: 0.126 seconds
%@  % MaxN = 4 ...   % CPU time: 0.214 seconds
%@  % MaxN = 5 ...   % CPU time: 0.331 seconds
%@  % MaxN = 6 ...   % CPU time: 0.455 seconds
%@  % MaxN = 7 ...   % CPU time: 0.634 seconds
%@  % MaxN = 8 ...   % CPU time: 0.808 seconds
%@  % MaxN = 9 ...   % CPU time: 1.009 seconds
%@  % MaxN = 10 ...   % CPU time: 1.268 seconds
%@  % MaxN = 11 ...   % CPU time: 1.514 seconds
%@  % MaxN = 12 ...   % CPU time: 1.818 seconds
%@ false.

%?- inconceivable((tally_pair(Q1, Q2, MaxN), qcompare(<, Q1, Q2), qcompare(~~, Q1, Q2)), MaxN, 1..12).
%@  % MaxN = 1 ...   % CPU time: 0.026 seconds
%@  % MaxN = 2 ...   % CPU time: 0.067 seconds
%@  % MaxN = 3 ...   % CPU time: 0.135 seconds
%@  % MaxN = 4 ...   % CPU time: 0.251 seconds
%@  % MaxN = 5 ...   % CPU time: 0.381 seconds
%@  % MaxN = 6 ...   % CPU time: 0.572 seconds
%@  % MaxN = 7 ...   % CPU time: 0.847 seconds
%@  % MaxN = 8 ...   % CPU time: 1.125 seconds
%@  % MaxN = 9 ...   % CPU time: 1.555 seconds
%@  % MaxN = 10 ...   % CPU time: 1.992 seconds
%@  % MaxN = 11 ...   % CPU time: 2.539 seconds
%@  % MaxN = 12 ...   % CPU time: 3.055 seconds
%@ false.

%% Show that =< and >= hold simultaneously only in case of equivalence:
%?- inconceivable((tally_pair(Q1, Q2, MaxN), qcompare(>=, Q1, Q2), qcompare(=<, Q1, Q2), dif(Q1, Q2)), MaxN, 1..12).
%@  % MaxN = 1 ...   % CPU time: 0.044 seconds
%@  % MaxN = 2 ...   % CPU time: 0.106 seconds
%@  % MaxN = 3 ...   % CPU time: 0.187 seconds
%@  % MaxN = 4 ...   % CPU time: 0.296 seconds
%@  % MaxN = 5 ...   % CPU time: 0.424 seconds
%@  % MaxN = 6 ...   % CPU time: 0.571 seconds
%@  % MaxN = 7 ...   % CPU time: 0.753 seconds
%@  % MaxN = 8 ...   % CPU time: 0.954 seconds
%@  % MaxN = 9 ...   % CPU time: 1.178 seconds
%@  % MaxN = 10 ...   % CPU time: 1.424 seconds
%@  % MaxN = 11 ...   % CPU time: 1.700 seconds
%@  % MaxN = 12 ...   % CPU time: 1.997 seconds
%@ false.

%% TODO: Refine the console output from inconceivable/3, with an understanding
%%       that it will generally be called on a Query that MUST FAIL.
%%       There may even be scope for omitting the separate Var argument
%%       by abstracting it automatically from the Query itself, recognized
%%       perhaps as the only unbound named variable.
inconceivable(Query, Var, Range) :-
    Var in Range,
    indomain(Var),
    format(" % MaxN = ~d ...", [Var]),
    time(call(Query)).

%?- inconceivable(violate_transitivity(=<, Q1, Q2, Q3, MaxN), MaxN, 3..5).
%@  % MaxN = 3 ...   % CPU time: 6.012 seconds
%@  % MaxN = 4 ...   % CPU time: 16.696 seconds
%@  % MaxN = 5 ...   % CPU time: 41.477 seconds
%@ false.

%% TODO: Generalize tally_pair/3, tally_triple/4 to LISTS of tallies?
tally_triple(T1/N1, T2/N2, T3/N3, Size) :-
    tally_pair(T1/N1, T2/N2, Size),
    N3 in 1..N2, indomain(N3),
    T3 in 0..N3.

%?- Size in 1..3, indomain(Size), tally_triple(Q1, Q2, Q3, Size).
%@    Size = 1, Q1 = _A/1, Q2 = _B/1, Q3 = _C/1, clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(_C in 0..1)
%@ ;  Size = 2, Q1 = _A/2, Q2 = _B/1, Q3 = _C/1, clpz:(_A in 0..2), clpz:(_B in 0..1), clpz:(_C in 0..1)
%@ ;  Size = 2, Q1 = _A/2, Q2 = _B/2, Q3 = _C/1, clpz:(_A in 0..2), clpz:(_B in 0..2), clpz:(_C in 0..1)
%@ ;  Size = 2, Q1 = _A/2, Q2 = _B/2, Q3 = _C/2, clpz:(_A in 0..2), clpz:(_B in 0..2), clpz:(_C in 0..2)
%@ ;  Size = 3, Q1 = _A/3, Q2 = _B/1, Q3 = _C/1, clpz:(_A in 0..3), clpz:(_B in 0..1), clpz:(_C in 0..1)
%@ ;  Size = 3, Q1 = _A/3, Q2 = _B/2, Q3 = _C/1, clpz:(_A in 0..3), clpz:(_B in 0..2), clpz:(_C in 0..1)
%@ ;  Size = 3, Q1 = _A/3, Q2 = _B/2, Q3 = _C/2, clpz:(_A in 0..3), clpz:(_B in 0..2), clpz:(_C in 0..2)
%@ ;  Size = 3, Q1 = _A/3, Q2 = _B/3, Q3 = _C/1, clpz:(_A in 0..3), clpz:(_B in 0..3), clpz:(_C in 0..1)
%@ ;  Size = 3, Q1 = _A/3, Q2 = _B/3, Q3 = _C/2, clpz:(_A in 0..3), clpz:(_B in 0..3), clpz:(_C in 0..2)
%@ ;  Size = 3, Q1 = _A/3, Q2 = _B/3, Q3 = _C/3, clpz:(_A in 0..3), clpz:(_B in 0..3), clpz:(_C in 0..3).

%% Demonstrate the TRANSITIVITY of (&=<) and (&>=)
violate_transitivity(C, Q1, Q2, Q3, Size) :-
    tally_triple(Q1, Q2, Q3, Size),
    qcompare(C, Q1, Q2),
    qcompare(C, Q2, Q3),
    %% At this point Q1 (C) Q2 (C) Q3 holds, and now we
    %% ask whether it's possible that Q1 (C) Q3 DOESN'T:
    qcompare(C, Q1, Q3, false). % reification to the rescue!

%?- inconceivable(violate_transitivity(>=, Q1, Q2, Q3, Size), Size, 1..12).
%@  % MaxN = 1 ...   % CPU time: 0.125 seconds
%@  % MaxN = 2 ...   % CPU time: 0.347 seconds
%@  % MaxN = 3 ...   % CPU time: 0.691 seconds
%@  % MaxN = 4 ...   % CPU time: 1.156 seconds
%@  % MaxN = 5 ...   % CPU time: 1.715 seconds
%@  % MaxN = 6 ...   % CPU time: 2.449 seconds
%@  % MaxN = 7 ...   % CPU time: 3.286 seconds
%@  % MaxN = 8 ...   % CPU time: 4.258 seconds
%@  % MaxN = 9 ...   % CPU time: 5.540 seconds
%@  % MaxN = 10 ...   % CPU time: 7.315 seconds
%@  % MaxN = 11 ...   % CPU time: 8.463 seconds
%@  % MaxN = 12 ...   % CPU time: 9.648 seconds
%@ false.

%?- inconceivable(violate_transitivity(=<, Q1, Q2, Q3, Size), Size, 1..12).
%@  % MaxN = 1 ...   % CPU time: 0.123 seconds
%@  % MaxN = 2 ...   % CPU time: 0.365 seconds
%@  % MaxN = 3 ...   % CPU time: 0.728 seconds
%@  % MaxN = 4 ...   % CPU time: 1.201 seconds
%@  % MaxN = 5 ...   % CPU time: 1.776 seconds
%@  % MaxN = 6 ...   % CPU time: 2.522 seconds
%@  % MaxN = 7 ...   % CPU time: 3.458 seconds
%@  % MaxN = 8 ...   % CPU time: 4.297 seconds
%@  % MaxN = 9 ...   % CPU time: 5.410 seconds
%@  % MaxN = 10 ...   % CPU time: 6.631 seconds
%@  % MaxN = 11 ...   % CPU time: 8.037 seconds
%@  % MaxN = 12 ...   % CPU time: 9.447 seconds
%@ false.

/*

Drawing inspiration from BOIN [2] as a specific example of cumulative-cohort
designs, we posit that in CCDs in general, dose-escalation decisions are
driven solely from BOUNDARY-HITTING EVENTS captured by =< and >= relations.
The boundaries are of 2 types: FLOORS that (when hit) indicate escalation,
and CEILINGS that indicate de-escalation or even dose removal when hit.

(I will prefer the terms FLOOR and CEILING because 'boundary' may suggest
a too-narrow 'topological' interpretation. We consider the boundary to
have been 'hit' (past participle!) even if the current tally lies INSIDE
the territory it bounds. If an accident at home results in a projectile
being embedded 1" into the ceiling or floor, rather than stuck at the
very surface, we still say the surface was hit.)

Arguably, sensible toxicity boundaries in dose-finding must take the form
of 'stairs' ascending from left to right in the T-N lattice depicted above.
Furthermore, adherence to 'reachability logic' effectively excludes step
heights greater than 1. (PROOF: Reachability effectively 'fills the corner'
in any proposed step of height 2, and this argument applies recursively to
higher steps.)

Thus all reasonable toxicity floors/ceilings can be constructed by a union
of 135-degree sectors generated by relations (_ &=< Q) and (_ &>= Q).
These unions are naturally represented as lists of the 'vertex tallies' Q.

Whereas step *heights* are always 1, the step *lengths* will generally be
longer---on the order of the reciprocal of the target toxicity rate (TTR).
Thus, we need not list a vertex tally for every denominator to define our
toxicity floors and ceilings.

For example, the decision boundaries from Liu & Yuan (2015) Table 1 ...

             cumulative patients treated at current dose (n_j):
              1    2    3    4    5    6    7    8    9    10    11    12

lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12

lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12

elimination   -    -   3/3  3/4  3/5  4/6  4/7  4/8  5/9  5/10  6/11  6/12

... may easily be represented by short lists:

  Lambda_1 = [0/1, 1/8]

  Lambda_2 = [1/1, 2/4, 3/6, 4/8, 5/11, 6/12]

  Eliminate = [3/5, 4/8, 5/10, 6/12]

This is because, e.g.,

  maplist((Q &=<), [0/1, 1/8]) <==> maplist((Q &=<), [0/1,...,0/7,1/8,...]),

which is a consequence of the TRANSITIVITY of (&=<):

  Q =< 0/1 ==> Q =< 0/N =< 0/1 for N > 1;
  Q =< 1/8 ==> Q =< 1/N =< 1/8 for N > 8.

*/

hit_ceiling_t(_, [], false).
hit_ceiling_t(Q, [C|Cs], Truth) :-
    if_(Q &>= C % this goal is deterministic for Q, C ∊ ℚ
	, Truth = true
	, hit_ceiling_t(Q, Cs, Truth)
       ).

hit_floor_t(_, [], false).
hit_floor_t(Q, [F|Fs], Truth) :-
    if_(Q &=< F % this goal is deterministic for Q, F ∊ ℚ
	, Truth = true
	, hit_floor_t(Q, Fs, Truth)
       ).


/*

Our 'spatial' understanding of floors/ceilings defined with respect to
(&=) and (&>=) shows that their minimal representation is obtained by
listing only their (outer) vertices.

*/

ceiling_vertex_t(Qs, T/N, Truth) :-
    member(T/N, Qs),
    N1 #= N + 1,
    if_(hit_ceiling_t(T/N1, Qs)
	, Truth = false
	, ( T_1 #= T - 1,
	    N_1 #= N - 1,
	    if_(hit_ceiling_t(T_1/N_1, Qs)
		, Truth = false
		, Truth = true
	       )
	  )
       ).

%?- ceiling_vertex_t([3/3,3/4,3/5,4/6,4/7,4/8,5/9,5/10,6/11,6/12], V, true).
%@    V = 3/5
%@ ;  V = 4/8
%@ ;  V = 5/10
%@ ;  V = 6/12
%@ ;  false.

%?- ceiling_vertex_t([3/3,3/4,3/5,4/6,5/7,6/8,6/11,6/12], V, true).
%@    V = 3/5
%@ ;  V = 6/12
%@ ;  false.

floor_vertex_t(Qs, T/N, Truth) :-
    member(T/N, Qs),
    N_1 #= N - 1,
    if_(hit_floor_t(T/N_1, Qs)
	, Truth = false
	, ( T1 #= T + 1,
	    N1 #= N + 1,
	    if_(hit_floor_t(T1/N1, Qs)
		, Truth = false
		, Truth = true
	       )
	  )
       ).

%?- floor_vertex_t([0/3,1/4,2/5,4/8,5/12], V, true).
%@    V = 2/5
%@ ;  V = 4/8
%@ ;  V = 5/12
%@ ;  false.

%?- floor_vertex_t([0/3,1/5,4/8,5/12], V, true).
%@    V = 0/3
%@ ;  V = 4/8
%@ ;  V = 5/12
%@ ;  false.

%?- floor_vertex_t([0/3,1/5,4/8,5/12], V, false).
%@    V = 1/5
%@ ;  false.

%?- tfilter(floor_vertex_t([0/3,1/5,4/8,5/12]), [0/3,1/5,4/8,5/12], Vs).
%@    Vs = [0/3,4/8,5/12]
%@ ;  false.

%% It will help to have unique, minimal ('canonical') representations
%% for floor- and ceiling-type boundaries.
%% TODO: Are these concepts better represented through SETS than lists?
%%       Shouldn't I require that each tally in a canonical floor/ceiling
%%       contributes something? Does this requirement ensure uniqueness?
ceiling_canonical(Qs, Ks) :-
    tfilter(ceiling_vertex_t(Qs), Qs, Ks_),
    sort(Ks_, Ks).

%?- ceiling_canonical([3/3,3/4,3/5,4/6,4/7,4/8,5/9,5/10,6/11,6/12], K).
%@    K = [3/5,4/8,5/10,6/12]
%@ ;  false.

%?- ceiling_canonical([6/12,4/8,3/3,3/4,3/5,4/6,4/7,5/9,5/10,6/11], K).
%@    K = [3/5,4/8,5/10,6/12]
%@ ;  false.

floor_canonical(Qs, Ks) :-
    tfilter(floor_vertex_t(Qs), Qs, Ks_),
    sort(Ks_, Ks).

%?- floor_canonical([0/1, 0/2, 0/3, 0/4, 0/5, 0/6, 0/7, 1/8, 1/9, 1/10, 1/11, 1/12], K).
%@    K = [0/1,1/8]
%@ ;  false.

%% tally_decision_ccd(?Q, ?Decision, +CCD) relates tallies Q to Decisions,
%% for a GIVEN ground cumulative-cohort design (CCD) which takes the form
%% of a triplet of boundaries, followed by max enrollment per cohort.
%% TODO: I believe the if_/3 cascade, with its default final 'escape clause',
%%       together with the determinism of hit_ceiling_t/3 and hit_floor_t/3,
%%       proves the determinism of tally_decision/2 so long as the defining
%%       boundaries are lists from ℚ.
tally_decision_ccd(Q, Decision, ccd(RemovalBdy, DeescBdy, EscBdy, FullCoh)) :-
    Q = T/N,
    N in 0..FullCoh, indomain(N),
    T in 0..N,
    if_(hit_ceiling_t(Q, RemovalBdy)
	, Decision = remove
	, if_(hit_ceiling_t(Q, DeescBdy)
	      , Decision = deescalate
	      , if_(hit_floor_t(Q, EscBdy)
		    , Decision = escalate
		    , Decision = stay
		   )
	     )
       ).

%% For testing purposes, we hard-code the CCD to obtain tally_decision/2
tally_decision(Q, Decision) :-
    tally_decision_ccd(Q, Decision, ccd([3/5, 4/8, 5/10, 6/12],
					[1/1, 2/4, 3/6, 4/8, 5/11],
					[0/1, 1/8],
					12)).

%?- tally_decision(T/5, Decision).
%@    Decision = stay, clpz:(T in 1..2)
%@ ;  Decision = escalate, T = 0
%@ ;  Decision = remove, clpz:(T in 3..5).

%?- tally_decision(Q, Decision).
%@    Q = 0/0, Decision = stay
%@ ;  Q = 0/1, Decision = escalate
%@ ;  Q = 1/1, Decision = deescalate
%@ ;  Q = 1/2, Decision = stay
%@ ;  Q = 0/2, Decision = escalate
%@ ;  Q = 2/2, Decision = deescalate
%@ ;  Q = 1/3, Decision = stay
%@ ;  Q = 0/3, Decision = escalate
%@ ;  Q = 2/3, Decision = deescalate
%@ ;  Q = 3/3, Decision = remove
%@ ;  Q = 1/4, Decision = stay
%@ ;  Q = 0/4, Decision = escalate
%@ ;  Q = 2/4, Decision = deescalate
%@ ;  Q = _A/4, Decision = remove, clpz:(_A in 3..4)
%@ ;  Q = _A/5, Decision = stay, clpz:(_A in 1..2)
%@ ;  Q = 0/5, Decision = escalate
%@ ;  ...


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Pick up from here...

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

/* An assumption we make about the dose-escalation trials modeled here
 * is that AT ANY GIVEN MOMENT they always partition the available doses
 * into 2 groups:
 * 'L' - A descending list of doses considered 'too Low' to enroll RIGHT NOW
 * 'R' - An ascending list of doses of which the head looks Right to enroll.
 * We can imagine these lists as 2 stacks, on the 'left' and 'right'.
 * 
 * I abandon the earlier 'bipartite' trial in alternating phases, in favor
 * of an all-in-one enroll/assess/decide cycle.
 *
 * it also appears necessary to append a 3rd, topmost list of excluded doses.
 */

/*
Above text is included for comparison with current outlook.
*/

%% Instead of carrying the tox boundaries along as parameters,
%% let's simply assert them into the database ...
%% Note that this naturally invites consideration of a DSL
%% in which such design rules could be expressed generally!

%% NB: The left-hand list of Ls ^ Rs is sorted in descending order.
%%     So heads L and R in [L|Ls] ^ [R|Rs] are tallies belonging to
%%     *adjacent* doses, notwithstanding their non-juxtaposition in
%%     our left-right reading of the term.

%% This enroll/3 goal, with its reified 'success' arg #3, creates fine
%% opportunities to bring various CCD-adapted STOPPING CRITERIA to bear.
%% Presently, we are simply checking whether cohort N0 is already 'full'.
enroll(T0/N0, T1/N1, Truth) :-
    N0 #>= 6 #<==> CohortFull, % suggestive of a line from a trial 'config file'?
    if_(CohortFull #= 1
	, Truth = false
	, ( N1 #= N0 + 1,
	    T in 0..1, % T is the pending tox assessment of newly-enrolled patient
	    indomain(T), % TODO: How to 'parametrize' this? Use OPTIONS?
	    T1 #= T0 + T,
	    Truth = true
	  )
       ).

%?- enroll(3/6, Q, Success).
%@    Q = _B/7, Success = true, clpz:(3+_A#=_B), clpz:(_A in 0..1), clpz:(_B in 3..4).

length_plus_1(Ls, MTD) :-
    length(Ls, MTD_1),
    MTD #= MTD_1 + 1.

stay(Ls ^ [R | Rs] ^ Es, State) :-
    if_(enroll(R, R1)
	, State = Ls ^ [R1 | Rs] ^ Es
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
	).

%% TODO: Simply deferring to the 'stay' case discards the information
%%       contained in having *desired* an escalation. Even if this
%%       accords with the 'standard' treatment of BOIN, it seems
%%       lacking in practical relevance. For all practical purposes,
%%       is it not the case that we would wish to impose a special
%%       escape clause so we don't have to run the top dose to 0/12
%%       (say) before declaring it as RP2D?
%% TODO: OTOH, does an analogy with Prolog program slicing hold lessons
%%       for us, about how to analyze the introduction of additional
%%       termination rules? Is there some sense in which such extra
%%       rules 'slice' a trial, where the paths are 'solutions'?
%%       Does this analogy advise AGAINST introducing an complicating
%%       logic like abovementioned "special escape clause"?
escalate(Ls ^ [R] ^ Es, State) :- % NB: this is a 'clamped' situation
    stay(Ls ^ [R] ^ Es, State).

escalate(Ls ^ [Q, R | Rs] ^ Es, State) :-
    if_(enroll(R, R1)
	, State = [Q | Ls] ^ [R1 | Rs] ^ Es
	%% If the next dose up (R) cannot be enrolled, that's because
	%% it's already full. What's more, it must have recommended
	%% de-escalation---which is how we got to the current dose Q!
	%% Accordingly, we declare the current dose to be MTD:
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
       ).

%% TODO: Suppose trial starts at lowest dose, and first patient has
%%       a DLT. If Deesc = [1/1,...] then the following rule would
%%       stop the trial right there! On the one hand, this does seem
%%       fully consistent with preference for safety, since initial
%%       dose has 'surprised' us. (So, we must rethink design, etc.)
%%       But OTOH this effectively turns an initial toxicity into a
%%       dose-removing event, blurring what otherwise ought to be a
%%       clear distinction between de-escalation and removal.
deescalate([] ^ _ ^ _, declare_mtd(0)). % deescalate from already-lowest dose

deescalate([L | Ls] ^ Rs ^ Es, State) :-
    if_(enroll(L, L1)
	, State = Ls ^ [L1 | Rs] ^ Es
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
       ).

remove(Ls ^ Rs ^ Es, State) :-
    append(Rs, Es, RsEs),
    deescalate(Ls ^ [] ^ RsEs, State).

%% TODO: As a purely practical matter, it may be necessary to implement
%%       STOPPING RULES that depend on TOTAL ENROLLMENT, and not merely
%%       upon CURRENT-COHORT enrollment. Otherwise, CPE may get bogged
%%       down in listing enormous numbers of highly unlikely paths.
%%       To preserve the pure CC-dependence of most of the predicates,
%        such a branch may have to be introduced here.
state0_action_state(Ls ^ [R | Rs] ^ Es, Action, State) :-
    %% TODO: Check whether total enrollment has been reached, and stop?
    tally_decision(R, Action), % NB: tally_decision/2 confers its determinism on the trial
    call(Action, Ls ^ [R | Rs] ^ Es, State).

%% Note how this state-machine naturally starts up from a blank slate:
%?- state0_action_state([] ^ [0/0, 0/0, 0/0] ^ [], Action, State).
%@    Action = stay, State = []^[_A/1,0/0,0/0]^[], clpz:(_A in 0..1). % now det! (vs below)
%@    Action = stay, State = []^[_A/1,0/0,0/0]^[], clpz:(_A in 0..1)
%@ ;  false.

%% TODO: Consider restoring an 'mtd_notfound' concept to the trial.
%%       Even if this notion disappears 'WLOG' from a purely formal perspective,
%%       it remains part of the lingo, and with good reason. There really IS a
%%       difference between having found a moderate number of DLTs at the RP2D,
%%       and having probed nowhere into the toxic dose region.
%% TODO: But consider whether this type of linguistic expansion trespasses into a
%%       strictly *pharmacologic* realm which a formal analysis might best avoid.
%%       A good test may be to demonstrate that mtd_notfound truly adds something
%%       essential that can't be 'tacked on' in a post-processing step.
%%actions(mtd_notfound(_)) --> [].
actions(declare_mtd(_)) --> [].
actions(S0) --> [(S0->A->S)],
		{ state0_action_state(S0, A, S) },
		actions(S).

%% Examine the smallest possible trial -- a trial with just 1 dose!
%?- phrase(actions([]^[0/0]^[]), Trial).
%@    Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[1/2]^[]),([]^[1/2]^[]->stay->[]^[1/3]^[]),([]^[1/3]^[]->stay->[]^[1/4]^[]),([]^[1/4]^[]->stay->[]^[_A/5]^[]),([]^[_A/5]^[]->stay->[]^[_C/6]^[]),([]^[_C/6]^[]->stay->[]^[... / ...]^[]),([]^[...]^[]->stay->[]^ ... ^[]),([]^ ... ^ ... ->stay-> ... ^ ...),(... -> ...)|...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(_E+_F#=_G), clpz:(_G+_H#=_I), clpz:(_I+_J#=_K), clpz:(_K+_L#=_M), clpz:(_M+_N#=_O), clpz:(1+_P#=_A), clpz:(_P in 0..1), clpz:(_A in 1..2), clpz:(_B in 0..1), clpz:(_C in 1..2), clpz:(_D in 0..1), clpz:(_E in 1..3), clpz:(_F in 0..1), clpz:(_G in 2..3), clpz:(_H in 0..1), clpz:(_I in 2..4), clpz:(_J in 0..1), clpz:(_K in 2..4), clpz:(_L in 0..1), clpz:(_M in 2..4), clpz:(_N in 0..1), clpz:(_O in 2..5)
%@ ;  Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[1/2]^[]),([]^[1/2]^[]->stay->[]^[1/3]^[]),([]^[1/3]^[]->stay->[]^[1/4]^[]),([]^[1/4]^[]->stay->[]^[_A/5]^[]),([]^[_A/5]^[]->stay->[]^[_C/6]^[]),([]^[_C/6]^[]->stay->[]^[... / ...]^[]),([]^[...]^[]->stay->[]^ ... ^[]),([]^ ... ^ ... ->stay-> ... ^ ...),(... -> ...)|...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(_E+_F#=_G), clpz:(_G+_H#=_I), clpz:(_I+_J#=4), clpz:(1+_K#=_A), clpz:(_K in 0..1), clpz:(_A in 1..2), clpz:(_B in 0..1), clpz:(_C in 1..2), clpz:(_D in 0..1), clpz:(_E in 1..3), clpz:(_F in 0..1), clpz:(_G in 2..3), clpz:(_H in 0..1), clpz:(_J in 0..1), clpz:(_I in 3..4)
%@ ;  Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[1/2]^[]),([]^[1/2]^[]->stay->[]^[1/3]^[]),([]^[1/3]^[]->stay->[]^[1/4]^[]),([]^[1/4]^[]->stay->[]^[_A/5]^[]),([]^[_A/5]^[]->stay->[]^[_C/6]^[]),([]^[_C/6]^[]->stay->[]^[... / ...]^[]),([]^[...]^[]->stay->[]^ ... ^[]),([]^ ... ^ ... ->stay-> ... ^ ...),(... -> ...)|...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(_E+_F#=3), clpz:(1+_G#=_A), clpz:(_G in 0..1), clpz:(_A in 1..2), clpz:(_B in 0..1), clpz:(_C in 1..2), clpz:(_D in 0..1), clpz:(_F in 0..1), clpz:(_E in 2..3)
%@ ;  Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[1/2]^[]),([]^[1/2]^[]->stay->[]^[1/3]^[]),([]^[1/3]^[]->stay->[]^[1/4]^[]),([]^[1/4]^[]->stay->[]^[1/5]^[]),([]^[1/5]^[]->stay->[]^[1/6]^[]),([]^[1/6]^[]->stay->[]^[... / ...]^[]),([]^[...]^[]->stay->[]^ ... ^[]),([]^ ... ^ ... ->escalate-> ... ^ ...),(... -> ...)|...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(2+_F#=_A), clpz:(_F in 0..1), clpz:(_A in 2..3), clpz:(_B in 0..1), clpz:(_C in 2..4), clpz:(_D in 0..1), clpz:(_E in 2..5)
%@ ;  ...

%% TODO: If this base case has (as I suspect) a 'closed-form' solution,
%%       then this opens up opportunities to try using Prolog to PROVE
%%       properties of cumulative-cohort designs.
%%       Might PROOFS BY INDUCTION be possible? Even if these proofs
%%       were of high time-complexity, the finite domain of the truly
%%       practical dose-escalation trials might even possess provable
%%       properties 'for all practical purposes'. This is the kind of
%%       result that invokes the spirit of Markus's exhaustive treatment
%%       of SGP.

%% TODO: Develop a 'condensed' format as done in 'aliquots.pl'
%% NB: Achieving an application of the previous condensed form
%%     (with any required extensions) will constitute progress
%%     toward comprehending 3+3 together with CCDs -- and maybe
%%     demonstrating their unification.
%% TODO: Suppose I collapse 'stay' actions?

/*

I wonder if CCDs have their own natural condensed representations.

First, consider how a pretty graphical layout would look.

For a trial that starts at dose 2:

0/0 0/1 0/0 0/0
     |
0/0 1/2 0/0 0/0
   /
0/1 1/2 0/0 0/0
 |
0/2 1/2 0/0 0/0
 |
0/3 1/2 0/0 0/0
   \
0/3 1/3 0/0 0/0

etc.

The diagram is clearer if we carry non-current tallies forward implicitly:

0/0 0/1 0/0 0/0
     |
    1/2
   /
0/1
 |
0/2
 |
0/3
   \
    1/3

Finally, the right or left margin could detail the binding constraint(s)
according to which each dose-escalation decision ( |, / or \) was made.
An attempt to escalate past top dose could be depicted with > in place of |.

-----

These visualizations also immediately suggest equivalent condensed string
expressions. Indeed, this can be done almost trivially by swapping / => :,
\ => ^, and | => -, in the implicit-LOCF form of the visualization.

For example, the above diagram immediately becomes:

(2) 0/1 - 1/2 : 0/1 - 0/2 - 0/3 ^ 1/3  

where the (2) indicates starting from dose-level 2. Alternatively, the
fully initial state of the trial could be specified instead, in an even
more direct transposition from the visualization format:

[0/0, 0/1, 0/0, 0/0] - 1/2 : 0/1 - 0/2 - 0/3 ^ 1/3.

-----

Other formats also suggest themselves as possibly more readable. Consider:

0/1@2 - 1/2@2 : 0/1@1 - 0/2@1 - 0/3@1 ^ 1/3@2,

or even

0/1@2 -x 1/2@2 :, 0/1@1 -o 0/2@1 -o 0/3@1 ^x 1/3@2,

which incorporates annotations for toxicity { o => 0, x => 1 }.

*/

%% Here's some code from aliquots.pl, to work from ...

%% Initially, I dispense with the toxicity indicators o & x,
%% since these add the complication of having to carry forward
%% the full set of tallies (at all doses) for differencing.

:- op(900, xfx, @).

condensed, [-, T/N@D] -->
    [ (_^_->stay->Ls^[(T/N)|_]^_) ],
    { length_plus_1(Ls, D) },
    condensed.
condensed, [^, T/N@D] -->
    [ (_^_->escalate->Ls^[(T/N)|_]^_) ],
    { length_plus_1(Ls, D) },
    condensed.
condensed, [:, T/N@D] -->
    [ (_^_->deescalate->Ls^[(T/N)|_]^_) ],
    { length_plus_1(Ls, D) },
    condensed.

condensed, [-, declare_mtd(MTD)] --> [ (_->stay->declare_mtd(MTD)) ].
condensed, [^, declare_mtd(MTD)] --> [ (_->escalate->declare_mtd(MTD)) ].
condensed, [v, declare_mtd(MTD)] --> [ (_->deescalate->declare_mtd(MTD)) ].
%%condensed --> []. %% Uncomment this for 'catch-all' permitting partial translations.

%% All the better to see you with ...
boin(Ls^Rs^Es, Translation) :- % BOIN trial starting from arbitrary state
    phrase(actions(Ls^Rs^Es), Trial),
    phrase(condensed, Trial, Translation).

boin(D, Path) :- % BOIN trial with D doses, starting at dose 1 by default
    length(Tallies, D),
    maplist(=(0/0), Tallies),
    boin([]^Tallies^[], Path).


%?- boin(1, Path), portray_clause(Path).
%@ [-,0/1@1,^,1/2@1,-,1/3@1,-,1/4@1,-,A/5@1,-,B/6@1,-,C/7@1,-,D/8@1,-,E/9@1,-,F/10@1,-,G/11@1,-,H/12@1,-,declare_mtd(1)].
%@    Path = [-,0/1@1,^,1/2@1,-,1/3@1,-,1/4@1,-,_A/5@1,-,_C/6@1,-,...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(_E+_F#=_G), clpz:(_G+_H#=_I), clpz:(_I+_J#=_K), clpz:(_K+_L#=_M), clpz:(_M+_N#=_O), clpz:(1+_P#=_A), clpz:(_P in 0..1), clpz:(_A in 1..2), clpz:(_B in 0..1), clpz:(_C in 1..2), clpz:(_D in 0..1), clpz:(_E in 1..3), clpz:(_F in 0..1), clpz:(_G in 2..3), clpz:(_H in 0..1), clpz:(_I in 2..4), clpz:(_J in 0..1), clpz:(_K in 2..4), clpz:(_L in 0..1), clpz:(_M in 2..4), clpz:(_N in 0..1), clpz:(_O in 2..5)
%@ ;  [-,0/1@1,^,1/2@1,-,1/3@1,-,1/4@1,-,A/5@1,-,B/6@1,-,C/7@1,-,D/8@1,-,E/9@1,-,4/10@1,-,5/11@1,v,declare_mtd(0)].
%@ Path = [-,0/1@1,^,1/2@1,-,1/3@1,-,1/4@1,-,_A/5@1,-,_C/6@1,-,...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(_E+_F#=_G), clpz:(_G+_H#=_I), clpz:(_I+_J#=4), clpz:(1+_K#=_A), clpz:(_K in 0..1), clpz:(_A in 1..2), clpz:(_B in 0..1), clpz:(_C in 1..2), clpz:(_D in 0..1), clpz:(_E in 1..3), clpz:(_F in 0..1), clpz:(_G in 2..3), clpz:(_H in 0..1), clpz:(_J in 0..1), clpz:(_I in 3..4)
%@ ;  ...

%?- boin(2, Path), portray_clause(Path).
%@ [-,0/1@1,^,0/1@2,^,1/2@2,-,1/3@2,-,1/4@2,-,A/5@2,-,B/6@2,-,C/7@2,-,D/8@2,-,E/9@2,-,F/10@2,-,G/11@2,-,H/12@2,-,declare_mtd(2)].
%@    Path = [-,0/1@1,^,0/1@2,^,1/2@2,-,1/3@2,-,1/4@2,-,_A/5@2,-,...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(_E+_F#=_G), clpz:(_G+_H#=_I), clpz:(_I+_J#=_K), clpz:(_K+_L#=_M), clpz:(_M+_N#=_O), clpz:(1+_P#=_A), clpz:(_P in 0..1), clpz:(_A in 1..2), clpz:(_B in 0..1), clpz:(_C in 1..2), clpz:(_D in 0..1), clpz:(_E in 1..3), clpz:(_F in 0..1), clpz:(_G in 2..3), clpz:(_H in 0..1), clpz:(_I in 2..4), clpz:(_J in 0..1), clpz:(_K in 2..4), clpz:(_L in 0..1), clpz:(_M in 2..4), clpz:(_N in 0..1), clpz:(_O in 2..5)
%@ ;  [-,0/1@1,^,0/1@2,^,1/2@2,-,1/3@2,-,1/4@2,-,A/5@2,-,B/6@2,-,C/7@2,-,D/8@2,-,E/9@2,-,4/10@2,-,5/11@2,:,1/2@1,-,1/3@1,-,1/4@1,-,F/5@1,-,G/6@1,-,H/7@1,-,I/8@1,-,J/9@1,-,K/10@1,-,L/11@1,-,M/12@1,-,declare_mtd(1)].
%@ Path = [-,0/1@1,^,0/1@2,^,1/2@2,-,1/3@2,-,1/4@2,-,_A/5@2,-,...], clpz:(_A+_B#=_C), clpz:(_C+_D#=_E), clpz:(_E+_F#=_G), clpz:(_G+_H#=_I), clpz:(_I+_J#=4), clpz:(_K+_L#=_M), clpz:(_M+_N#=_O), clpz:(_O+_P#=_Q), clpz:(_Q+_R#=_S), clpz:(_S+_T#=_U), clpz:(_U+_V#=_W), clpz:(_W+_X#=_Y), clpz:(1+_Z#=_A), clpz:(1+_A1#=_K), clpz:(_Z in 0..1), clpz:(_A in 1..2), clpz:(_B in 0..1), clpz:(_C in 1..2), clpz:(_D in 0..1), clpz:(_E in 1..3), clpz:(_F in 0..1), clpz:(_G in 2..3), clpz:(_H in 0..1), clpz:(_J in 0..1), clpz:(_I in 3..4), clpz:(_A1 in 0..1), clpz:(_K in 1..2), clpz:(_L in 0..1), clpz:(_M in 1..2), clpz:(_N in 0..1), clpz:(_O in 1..3), clpz:(_P in 0..1), clpz:(_Q in 2..3), clpz:(_R in 0..1), clpz:(_S in 2..4), clpz:(_T in 0..1), clpz:(_U in 2..4), clpz:(_V in 0..1), clpz:(_W in 2..4), clpz:(_X in 0..1), clpz:(_Y in 2..5)
%@ ;  ...

%% I'd like to complete the translation to the T[,,] arrays
%% that support the formalism of WWTT <arXiv:2012.05301>.
/*

How did I do this for 3+3, and what (if anything) has to change?


The indexing of cohorts on a per-dose basis nicely carries over into CCDs.
So the individual matrices T[,,j], j in 1..J retain the same form and semantics.
Furthermore, with a deterministic trial design, the path could be 'read out'
from these matrices quite straightforwardly.

But our plan for a more intimate involvement of Prolog in package 'precautionary'
happily liberates us from the supposed need to explicitly represent every step of
trial operation to R. When we need to determine properties of the trial paths at
such levels of detail, clearly the appropriate setting for the investigations is
Prolog itself. Thus, we need not even preserve all of the path information in the
cached arrays T[,,].

For CCDs with cohorts of 1, all R really needs to obtain from these arrays are
dose-wise tallies of x's and o's. This information suffices for computing path
probabilities and ordinalized toxicity rates, and indeed even for extracting the
dose recommendations provided that the stopping rule is CC-adapted and known via
context such as indexes into the data structure storing the designs' T[,,] arrays.

Thus, all we really need to extract for T[,,] output of each design is the final
state of actions/1.

*/

:- op(900, xfx, ~>).

%% First, extract the path matrix:
path_matrix, [S ~> MTD] --> [ (S->_->declare_mtd(MTD)) ].
path_matrix, [] --> [(_->_->_^_^_)], path_matrix. % 'skip to end'

ccd_matrix(D, Matrix) :-
    length(Tallies, D), maplist(=(0/0), Tallies),
    phrase(actions([]^Tallies^[]), Path),
    phrase(path_matrix, Path, Matrix).

%% TODO: Can I label 'on demand', e.g. so that it happens right before portray_clause/2?
%%       Perhaps my solution is to redefine enroll/3 for each new design? If it turns out
%%       that everything variable in CCDs is inside enroll/3, that's a splendid finding,
%%       since it focuses attention for development of a CCD DSL.
%?- ccd_matrix(2, Matrix).
%@    Matrix = [[0/1]^[0/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/2,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/4,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/3,3/6]^[]~>0]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/2,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/4,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/3,3/6]^[]~>0]
%@ ;  Matrix = [[0/2]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/3,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/4,3/6]^[]~>0]
%@ ;  Matrix = [[]^[0/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/5]~>0]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/5]~>0]
%@ ;  Matrix = [[]^[2/4]^[3/5]~>0]
%@ ;  Matrix = [[]^[1/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/4]^[]~>0]
%@ ;  Matrix = [[]^[2/4,2/4]^[]~>0]
%@ ;  Matrix = [[]^[2/3,2/4]^[]~>0]
%@ ;  Matrix = [[0/1]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  Matrix = [[0/1]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/2,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/4,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/3,3/6]^[]~>0]
%@ ;  Matrix = [[0/2]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/3,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/4,3/6]^[]~>0]
%@ ;  Matrix = [[]^[0/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/5]~>0]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/5]~>0]
%@ ;  Matrix = [[]^[2/4]^[3/5]~>0]
%@ ;  Matrix = [[]^[1/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/4]^[]~>0]
%@ ;  Matrix = [[]^[2/4,2/4]^[]~>0]
%@ ;  Matrix = [[]^[2/3,2/4]^[]~>0]
%@ ;  Matrix = [[0/3]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/4,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[0/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/5]~>0]
%@ ;  Matrix = [[]^[1/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/4]^[]~>0]
%@ ;  Matrix = [[]^[2/4,2/4]^[]~>0]
%@ ;  Matrix = [[]^[0/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/4]~>0]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/4]~>0]
%@ ;  Matrix = [[]^[2/4]^[3/4]~>0]
%@ ;  Matrix = [[]^[1/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/3]^[]~>0]
%@ ;  Matrix = [[]^[2/4,2/3]^[]~>0]
%@ ;  Matrix = [[]^[2/3,2/3]^[]~>0]
%@ ;  Matrix = [[0/2]^[1/6]^[]~>2]
%@ ;  Matrix = [[0/2]^[2/6]^[]~>2]
%@ ;  Matrix = [[0/2]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/3,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[2/4,3/6]^[]~>0]
%@ ;  Matrix = [[0/3]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/4,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[3/6,3/6]^[]~>0]
%@ ;  Matrix = [[]^[0/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/5]~>0]
%@ ;  Matrix = [[]^[1/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/4]^[]~>0]
%@ ;  Matrix = [[]^[2/4,2/4]^[]~>0]
%@ ;  Matrix = [[0/4]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/5,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[2/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[0/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/4]^[]~>0]
%@ ;  Matrix = [[]^[0/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/4]~>0]
%@ ;  Matrix = [[]^[1/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/3]^[]~>0]
%@ ;  Matrix = [[]^[2/4,2/3]^[]~>0]
%@ ;  Matrix = [[0/5]^[2/6]^[]~>2]
%@ ;  Matrix = [[]^[0/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[1/6,3/6]^[]~>1]
%@ ;  Matrix = [[]^[0/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/5]~>1]
%@ ;  Matrix = [[]^[1/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/4]^[]~>1]
%@ ;  Matrix = [[]^[0/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/4]~>1]
%@ ;  Matrix = [[]^[1/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/3]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/3]^[]~>0]
%@ ;  Matrix = [[]^[0/6]^[3/3]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/3]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/3]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/3]~>1]
%@ ;  Matrix = [[]^[1/6]^[3/3]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/3]~>1]
%@ ;  Matrix = [[]^[2/6]^[3/3]~>1]
%@ ;  Matrix = [[]^[3/6]^[3/3]~>0]
%@ ;  Matrix = [[]^[1/6,2/2]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/2]^[]~>1]
%@ ;  Matrix = [[]^[2/6,2/2]^[]~>1]
%@ ;  Matrix = [[]^[3/6,2/2]^[]~>0]
%@ ;  Matrix = [[]^[2/4,2/2]^[]~>0]
%@ ;  Matrix = [[]^[1/6,1/1]^[]~>1]
%@ ;  Matrix = [[]^[2/6,1/1]^[]~>1]
%@ ;  Matrix = [[]^[2/6,1/1]^[]~>1]
%@ ;  Matrix = [[]^[3/6,1/1]^[]~>0]
%@ ;  Matrix = [[]^[2/4,1/1]^[]~>0]
%@ ;  Matrix = [[]^[2/3,1/1]^[]~>0]
%@ ;  Matrix = [[]^[1/1,0/0]^[]~>0]
%@ ;  false.

%?- findall(Matrix, ccd_matrix(2, Matrix), Paths), length(Paths, J).
%@    Paths = [[[0/1]^[0/6]^[]~>2],[[0/1]^[1/6]^[]~>2],[[0/1]^[1/6]^[]~>2],[[0/1]^[2/6]^[]~>2],[[0/1]^[1/6]^[]~>2],[[0/1]^[2/6]^[]~>2],[[... / ...]^[...]^[]~>2],[[]^ ... ^ ... ~>1],[... ~> ...],...|...], J = 212.

%?- use_module(library(lambda)).
%@    true.
%?- J+\(time(findall(Matrix, ccd_matrix(2, Matrix), _Paths)), length(_Paths, J)).
%@    % CPU time: 114.852 seconds
%@    % CPU time: 114.856 seconds
%@    J = 212.

%?- J+\(time(findall(Matrix, ccd_matrix(3, Matrix), _Paths)), length(_Paths, J)).
%@    % CPU time: 629.949 seconds
%@    % CPU time: 629.953 seconds
%@    J = 1151.

%?- J+\(time(findall(Matrix, ccd_matrix(4, Matrix), _Paths)), length(_Paths, J)).
%@    % CPU time: 3509.381 seconds
%@    % CPU time: 3509.385 seconds
%@    J = 6718.
