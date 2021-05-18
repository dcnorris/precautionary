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
:- use_module(library(time)).
:- use_module(library(debug)).
:- use_module(library(dif)).
:- use_module(library(reif)).
:- use_module(library(format)).
:- use_module(library(between)).

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

% My initial emphasis is on generating all possible paths (CPE) for the BOIN
% design set forth in the table above. Although the BOIN design of [1] lacks
% any terminating principle except elimination of all doses, we do need such
% a rule here. The most natural type of rule, in view of the Table above,
% might be a 'stop-for-consensus' type of rule as found in package 'dtpcrm'.
% This is specified as a maximum number of patients to enroll at any 1 dose.

%% In general, I may require a max-enrollment *parameter*,
%% however unsightly it may be tagging along like this ...
tally(DLTs/Enrolled, MaxN) :-
    ground(MaxN),
    Enrolled in 0..MaxN,
    indomain(Enrolled),
    DLTs in 0..Enrolled,
    indomain(DLTs).

%% Having a 'default' max cohort size will make dev & test a bit easier:
tally(DLTs/Enrolled) :- tally(DLTs/Enrolled, 12).

%% tally/1 terminates:
%?- tally(_), false.
%@ false.

%% -----------------------------------------------------------------------------

%% Is discontiuous directive not working?
:- discontiguous(qcompare/3). % [this is line 144 reference below]
%% caught: error(syntax_error(incomplete_reduction),read_term/3:144)
%% false.

%% The most fundamental relation between tallies is REACHABILITY,
%% the question of whether a possible path exists connecting them.
qcompare(~~, T1/N1, T2/N2) :- % REACHABILITY
    %% The following conditions translate as follows, under 3 scenarios:
    %% (=) If N1==N2, then we get T2 =< T1 =< T2 which holds iff T1==T2.
    %% (<) If N1 < N2, these conds translate to T1 =< T2 =< T1 + MaxTox.
    %% (>) If N1 > N2, these conds translate to T2 =< T1 =< T2 + MaxTox.
    T2 #=< T1 + max(0, N2 - N1),
    T1 #=< T2 + max(0, N1 - N2).

%% A helper to halve the results from the little query below...
earlier(_/N1, _/N2) :- N1 #< N2.

%% A 'general-enough' query, if not the MGQ:
%?- tally(Q1, 3), tally(Q2, 3), earlier(Q1, Q2), qcompare(~~, Q1, Q2).
%@    Q1 = 0/0, Q2 = 0/1
%@ ;  Q1 = 0/0, Q2 = 1/1
%@ ;  Q1 = 0/0, Q2 = 0/2
%@ ;  Q1 = 0/0, Q2 = 1/2
%@ ;  Q1 = 0/0, Q2 = 2/2
%@ ;  Q1 = 0/0, Q2 = 0/3
%@ ;  Q1 = 0/0, Q2 = 1/3
%@ ;  Q1 = 0/0, Q2 = 2/3
%@ ;  Q1 = 0/0, Q2 = 3/3
%@ ;  Q1 = 0/1, Q2 = 0/2
%@ ;  Q1 = 0/1, Q2 = 1/2
%@ ;  Q1 = 0/1, Q2 = 0/3
%@ ;  Q1 = 0/1, Q2 = 1/3
%@ ;  Q1 = 0/1, Q2 = 2/3
%@ ;  Q1 = 1/1, Q2 = 1/2
%@ ;  Q1 = 1/1, Q2 = 2/2
%@ ;  Q1 = 1/1, Q2 = 1/3
%@ ;  Q1 = 1/1, Q2 = 2/3
%@ ;  Q1 = 1/1, Q2 = 3/3
%@ ;  Q1 = 0/2, Q2 = 0/3
%@ ;  Q1 = 0/2, Q2 = 1/3
%@ ;  Q1 = 1/2, Q2 = 1/3
%@ ;  Q1 = 1/2, Q2 = 2/3
%@ ;  Q1 = 2/2, Q2 = 2/3
%@ ;  Q1 = 2/2, Q2 = 3/3
%@ ;  false.

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

%% Operators for the truly useful comparisons of BOIN, restricted to
%% the simplex of valid tallies:
:- op(900, xfx, &=<).
&=<(Q1, Q2) :-
    tally(Q1), tally(Q2),
    qcompare(=<, Q1, Q2).

&=<(Q1, Q2, Truth) :- % reified
    tally(Q1), tally(Q2),
    qcompare(=<, Q1, Q2, Truth).

:- op(900, xfx, &>=).
&>=(Q1, Q2) :-
    tally(Q1), tally(Q2),
    qcompare(>=, Q1, Q2).

&>=(Q1, Q2, Truth) :- % reified
    tally(Q1), tally(Q2),
    qcompare(>=, Q1, Q2, Truth).


%% Demonstrate that strict inequalities are exclusive of ~~
%?- inconceivable((tally(Q1, MaxN), tally(Q2, MaxN), qcompare(>, Q1, Q2), qcompare(~~, Q1, Q2)), MaxN, 6..12).
%@  % MaxN = 6 ...   % CPU time: 2.435 seconds
%@  % MaxN = 7 ...   % CPU time: 3.814 seconds
%@  % MaxN = 8 ...   % CPU time: 5.734 seconds
%@  % MaxN = 9 ...   % CPU time: 8.357 seconds
%@  % MaxN = 10 ...   % CPU time: 11.836 seconds
%@  % MaxN = 11 ...   % CPU time: 16.735 seconds
%@  % MaxN = 12 ...   % CPU time: 23.156 seconds
%@ false.

%?- inconceivable((tally(Q1, MaxN), tally(Q2, MaxN), qcompare(<, Q1, Q2), qcompare(~~, Q1, Q2)), MaxN, 1..12).
%@  % MaxN = 1 ...   % CPU time: 0.069 seconds
%@  % MaxN = 2 ...   % CPU time: 0.174 seconds
%@  % MaxN = 3 ...   % CPU time: 0.392 seconds
%@  % MaxN = 4 ...   % CPU time: 0.783 seconds
%@  % MaxN = 5 ...   % CPU time: 1.418 seconds
%@  % MaxN = 6 ...   % CPU time: 2.413 seconds
%@  % MaxN = 7 ...   % CPU time: 3.822 seconds
%@  % MaxN = 8 ...   % CPU time: 5.781 seconds
%@  % MaxN = 9 ...   % CPU time: 8.410 seconds
%@  % MaxN = 10 ...   % CPU time: 11.865 seconds
%@  % MaxN = 11 ...   % CPU time: 16.358 seconds
%@  % MaxN = 12 ...   % CPU time: 21.935 seconds
%@ false.

%% Show that =< and >= hold simultaneously only in case of equivalence:
%?- time((tally(Q1), tally(Q2), qcompare(>=, Q1, Q2), qcompare(=<, Q1, Q2), dif(Q1, Q2))).
%@    % CPU time: 22.008 seconds
%@ false.
%?- time((Q1 &>= Q2, Q1 &=< Q2, dif(Q1, Q2))).
%@    % CPU time: 35.625 seconds % TODO: Slower because (&=<) re-posts tally/1 constraints?
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

%% Demonstrate the TRANSITIVITY of (&=<) and (&>=)
violate_transitivity(C, Q1, Q2, Q3, MaxN) :-
    tally(Q1, MaxN),
    tally(Q2, MaxN),
    qcompare(C, Q1, Q2),
    tally(Q3, MaxN),
    qcompare(C, Q2, Q3),
    %% At this point Q1 (C) Q2 (C) Q3 holds, and now we
    %% ask whether it's possible that Q1 (C) Q3 DOESN'T:
    qcompare(C, Q1, Q3, false). % reification to the rescue!

%?- inconceivable(violate_transitivity(>=, Q1, Q2, Q3, MaxN), MaxN, 5..12).
%@  % MaxN = 5 ...   % CPU time: 40.463 seconds
%@  % MaxN = 6 ...   % CPU time: 87.731 seconds
%@  % MaxN = 7 ...   % CPU time: 171.186 seconds
%@  % MaxN = 8 ...   % CPU time: 321.167 seconds
%@  % MaxN = 9 ...   % CPU time: 553.290 seconds
%@  % MaxN = 10 ...   % CPU time: 919.145 seconds
%@  % MaxN = 11 ...   % CPU time: 1459.348 seconds
%@  % MaxN = 12 ...   % CPU time: 2246.696 seconds
%@ false.

%?- inconceivable(violate_transitivity(=<, Q1, Q2, Q3, MaxN), MaxN, 3..5).
%@  % MaxN = 3 ...   % CPU time: 6.333 seconds
%@  % MaxN = 4 ...   % CPU time: 18.381 seconds
%@  % MaxN = 5 ...   % CPU time: 42.346 seconds
%@ false.

/*

In BOIN (and perhaps CCDs generally?), dose-escalation decisions can be
driven solely from BOUNDARY-HITTING EVENTS captured by =< and >= relations.
The boundaries are of 2 types: FLOORS that (when hit) indicate escalation,
and CEILINGS that indicate de-escalation or even dose removal when hit.

(I will prefer the terms FLOOR and CEILING because 'boundary' may suggest
a too-narrow 'topological' interpretation. We consider the boundary to
have been 'hit' (past participle!) even if the current tally lies INSIDE
the territory it bounds. If an accident at home results in a projectile
being embedded 1" into the ceiling or floor, rather than stuck at the
very surface, we still say the surface was hit.)

Arguably, all sensible toxicity boundaries dose-finding will take the form
of 'stairs' ascending from left to right in the T-N lattice depicted above.
Furthermore, adherence to 'reachability logic' effectively excludes step
heights greater than 1. (TODO: Prove this.)

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
    if_(Q &>= C
	, Truth = true
	, hit_ceiling_t(Q, Cs, Truth)
       ).

hit_floor_t(_, [], false).
hit_floor_t(Q, [F|Fs], Truth) :-
    if_(Q &=< F
	, Truth = true
	, hit_floor_t(Q, Fs, Truth)
       ).


tally_decision(Q, Decision) :-
    RemovalBdy = [3/5, 4/8, 5/10, 6/12],   % To begin, we hard-code the
    DeescBdy = [1/1, 2/4, 3/6, 4/8, 5/11], % trial definition, according
    EscBdy = [0/1, 1/8],                   % to Table 1 abstracted above.
    tally(Q),
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

%?- tally_decision(2/5, Decision).
%@    Decision = stay.

%?- tally_decision(3/5, Decision).
%@    Decision = remove.

%?- tally_decision(1/5, Decision).
%@    Decision = stay.

%?- tally_decision(0/5, Decision).
%@    Decision = escalate.

%?- tally_decision(T/5, Decision).
%@    Decision = escalate, T = 0
%@ ;  Decision = stay, T = 1
%@ ;  Decision = stay, T = 2
%@ ;  Decision = remove, T = 3
%@ ;  Decision = remove, T = 4
%@ ;  Decision = remove, T = 5.

%?- tally_decision(Q, Decision).
%@    Q = 0/0, Decision = stay % <-- this one is interesting! TODO: consider
%@ ;  Q = 0/1, Decision = escalate
%@ ;  Q = 1/1, Decision = deescalate
%@ ;  Q = 0/2, Decision = escalate
%@ ;  Q = 1/2, Decision = stay
%@ ;  Q = 2/2, Decision = deescalate
%@ ;  Q = 0/3, Decision = escalate
%@ ;  Q = 1/3, Decision = stay
%@ ;  Q = 2/3, Decision = deescalate
%@ ;  Q = 3/3, Decision = remove
%@ ;  Q = 0/4, Decision = escalate
%@ ;  Q = 1/4, Decision = stay
%@ ;  Q = 2/4, Decision = deescalate
%@ ;  Q = 3/4, Decision = remove
%@ ;  Q = 4/4, Decision = remove
%@ ;  Q = 0/5, Decision = escalate
%@ ;  Q = 1/5, Decision = stay
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
 * Since it proves useful to describe the trial in two alternating phases,
 * we designate the state of the trial via L and R lists conjoined with
 * alternating infix functors:
 * '^' - Imagine: 'scales of Justice' where we make judgements about tallies
 * ':' - Imagine: an ellipsis (or dripping faucet) where we enroll patients
 *       (who 'trickle in'), wait for & assess toxicities.
 */

/*
Above text is included for comparison with current outlook.
*/

%% Instead of carrying the tox boundaries along as parameters,
%% let's simply assert them into the database ...
%% Note that this naturally invites consideration of a DSL
%% in which such design rules could be expressed generally!

escalate(Q) :- tally_decision(Q, escalate).
deescalate(Q) :- tally_decision(Q, deescalate).
remove(Q) :- tally_decision(Q, remove).
stay(Q) :- tally_decision(Q, stay).

%% NB: The left-hand list of Ls ^ Rs is sorted in descending order.
%%     So heads L and R in [L|Ls] ^ [R|Rs] are tallies belonging to
%%     *adjacent* doses, notwithstanding their non-juxtaposition in
%%     our left-right reading of the term.

%% The generic stopping condition obtains when an attempt is made
%% to enroll a 'full' cohort, at which time the trial stops with
%% that dose recommended as RP2D. (Even the seemingly 'special' case
%% of attempted deescalation from the lowest dose fits this pattern,
%% provided that we treat the zero dose as 'already full' -- which
%% in a pharmacologic sense it really is, since we know everything
%% about its effects already.)
cohort_full(_/N, TF) :- cohort_full(N, TF).
cohort_full(N, true) :- N #>= 12.
cohort_full(N, false) :- N in 0..11.

enroll(T0/N0, T1/N1, true) :-
    cohort_full(N0, false),
    N1 #= N0 + 1,
    T in 0..1, % T is the pending tox assessment of newly-enrolled patient
    T1 #= T0 + T.

enroll(Q0, _, false) :-
    cohort_full(Q0, true).

enroll(Q1, Q1) :- enroll(Q0, Q1, true).

%% Overload enroll/2 on lists of tallies, enrolling head tally
enroll([T0/N0 | Qs], [T1/N1 | Qs]) :-
    enroll(T0/N0, T1/N1).

length_plus_1(Ls, MTD) :-
    length(Ls, MTD_1),
    MTD #= MTD_1 + 1.

stay(Ls ^ [R | Rs], State) :-
    if_(enroll(R, R1)
	, State = Ls ^ [R1 | Rs]
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
	).
		
escalate(Ls ^ [R], State) :- % NB: this is a 'clamped' situation
    stay(Ls ^ [R], State).

escalate(Ls ^ [Q, R | Rs], State) :-
    if_(enroll(R, R1)
	, State = [Q | Ls] ^ [R1 | Rs]
	, ( length([R,Q|Ls], MTD),
	    declare_mtd(MTD)
	  )
       ).

deescalate([] ^ _, declare_mtd(0)). % deescalate from already-lowest dose

deescalate([L | Ls] ^ Rs, Ls ^ [L1 | Rs]) :- enroll(L, L1).
deescalate([L | Ls] ^ _, declare_mtd(MTD)) :- cohort_full(L, true),
					      length_plus_1(Ls, MTD).

remove(Ls ^ [R | Rs], State) :- deescalate(Ls ^ [], State).

%% TODO: Note how TRIVIAL state0_action_state/3 has become.
%%       Does this demand refactoring, or does it represent other
%%       types of opportunity -- even meta-interpretation?

state0_action_state(Ls ^ [R | Rs], escalate, State) :-
    escalate(R),
    escalate(Ls ^ [R | Rs], State).

state0_action_state(Ls ^ [R | Rs], deescalate, State) :-
    deescalate(R),
    deescalate(Ls ^ [R | Rs], State).

state0_action_state(Ls ^ [R | Rs], stay, State) :-
    stay(R),
    stay(Ls ^ [R | Rs], State).

state0_action_state(Ls ^ [R | Rs], remove, State) :-
    remove(R),
    remove(Ls ^ [R | Rs], State).

%% Note how this state-machine naturally starts up from a blank slate:
%?- state0_action_state([] ^ [0/0, 0/0, 0/0], Action, State).
%@    Action = stay, State = []^[_A/1,0/0,0/0], clpz:(_A in 0..1)
%@ ;  false.
%% NB: Without loading library(between), I obtain the following error:
%@ caught: error(existence_error(procedure,between/3),between/3)

%% TODO: Consider restoring an 'mtd_notfound' concept to the trial.
%%       Even if this notion disappears 'WLOG' from a purely formal perspective,
%%       it remains part of the lingo, and with good reason. There really IS a
%%       difference between having found a moderate number of DLTs at the RP2D,
%%       and having probed nowhere into the toxic dose region.
%%actions(mtd_notfound(_)) --> [].
actions(declare_mtd(_)) --> [].
actions(S0) --> [(S0->A->S)],
		{ state0_action_state(S0, A, S) },
		actions(S).

%% TODO: Try portray_clause/1 to visualize Trial below,
%%       in advance of a suitable condensed reformatting.

%% Examine the smallest possible trial -- a trial with just 1 dose!
%?- phrase(actions([]^[0/0]), Trial), portray_clause(Trial).
%@ [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[0/7]->escalate->[]^[0/8]),([]^[0/8]->escalate->[]^[0/9]),([]^[0/9]->escalate->[]^[0/10]),([]^[0/10]->escalate->[]^[0/11]),([]^[0/11]->escalate->[]^[0/12]),([]^[0/12]->escalate->declare_mtd(1))].
%@    Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[0/7]->escalate->[]^[0/8]),([]^[0/8]->escalate->[]^[0/9]),([]^[0/9]->escalate->[]^[0/10]),([]^[0/10]->escalate->[]^[0/11]),([]^[0/11]->escalate->[]^[1/12]),([]^[1/12]->escalate->declare_mtd(1))].
%@ Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[0/7]->escalate->[]^[0/8]),([]^[0/8]->escalate->[]^[0/9]),([]^[0/9]->escalate->[]^[0/10]),([]^[0/10]->escalate->[]^[1/11]),([]^[1/11]->escalate->[]^[1/12]),([]^[1/12]->escalate->declare_mtd(1))].
%@ Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
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
