% Abstracting the cumulative-cohort design (CCD) principle for dose escalation

:- module(ccd, [
	      ceiling_canonical/2,
	      floor_canonical/2,
	      tally_decision_ccd/3,
	      %% This infix op is used by path_matrix//0 to represent final rec:
              op(900, xfx, ~>),
	      cohort_max/1, % dynamic predicates that
	      enroll_max/1, % client code should redefine
	      ccd_d_matrix/3,
	      ccd_d_pathvector/3,
	      ccd_d_cmax_nmax_tabfile/5
	  ]).

:- dynamic(cohort_max/1).
:- dynamic(enroll_max/1).

:- multifile(cohort_max/1).
:- multifile(enroll_max/1).

/*

The 'cumulative-cohort design' (CCD) concept in [1] establishes a class of
dose-escalation designs in which escalation decisions are a function strictly
of the cumulative toxicity rate AT THE CURRENT DOSE. These designs are thus
positioned midway between what [1] calls 'group designs' (in which escalation
depends only on the rate observed in the current COHORT), and general designs
where these decisions depend on the observed rates AT ALL DOSES. (Note that,
in the dose-escalation context where patients are treated as *exchangeable*,
the list of dose-wise tallies [T1/N1, T2/N2, ..., Td/Nd] is a sufficient
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
:- use_module(tally).
:- use_module(benchmarking).

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
    if_(qcompare(>=, Q, C)
	, Truth = true
	, hit_ceiling_t(Q, Cs, Truth)
       ).

hit_floor_t(_, [], false).
hit_floor_t(Q, [F|Fs], Truth) :-
    if_(qcompare(=<, Q, F)
	, Truth = true
	, hit_floor_t(Q, Fs, Truth)
       ).


/*

Our 'spatial' understanding of floors/ceilings defined with respect to
(&=) and (&>=) shows that their minimal representation is obtained by
listing only their (outer) vertices.

*/

ceiling_vertex_t(Qs, T/N, Truth) :-
    memberd_t(T/N, Qs, true),
    N1 #= N + 1,
    if_(hit_ceiling_t(T/N1, Qs)
	, Truth = false
	, ( Tminus1 #= T - 1,
	    Nminus1 #= N - 1,
	    if_(hit_ceiling_t(Tminus1/Nminus1, Qs)
		, Truth = false
		, Truth = true
	       )
	  )
       ).

%?- time(ceiling_vertex_t([3/3,3/4,3/5,4/6,4/7,4/8,5/9,5/10,6/11,6/12], V, true)).
%@    % CPU time: 0.017 seconds
%@    V = 3/5
%@ ;  % CPU time: 0.061 seconds
%@    V = 4/8
%@ ;  % CPU time: 0.105 seconds
%@    V = 5/10
%@ ;  % CPU time: 0.154 seconds
%@    V = 6/12
%@ ;  % CPU time: 0.174 seconds
%@    false. %    ^^^^^ using memberd_t(T/N, Qs, true)
%@    % CPU time: 0.011 seconds
%@    V = 3/5
%@ ;  % CPU time: 0.031 seconds
%@    V = 4/8
%@ ;  % CPU time: 0.052 seconds
%@    V = 5/10
%@ ;  % CPU time: 0.074 seconds
%@    V = 6/12
%@ ;  % CPU time: 0.086 seconds
%@    false. %    ^^^^^ using member(T/N, Qs)

%?- ceiling_vertex_t([3/3,3/4,3/5,4/6,5/7,6/8,6/11,6/12], V, true).
%@    V = 3/5
%@ ;  V = 6/12
%@ ;  false.

floor_vertex_t(Qs, T/N, Truth) :-
    memberd_t(T/N, Qs, true),
    Nminus1 #= N - 1,
    if_(hit_floor_t(T/Nminus1, Qs)
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

cohort_full(N, Truth) :-
    cohort_max(Nmax),
    if_(N #< Nmax
	, Truth = false
	, Truth = true
       ).

enroll(T0/N0, T1/N1, Truth) :-
    if_(cohort_full(N0)
       , Truth = false
       , ( N1 #= N0 + 1,
           T in 0..1, % T is the pending tox assessment of newly-enrolled patient
           indomain(T), % TODO: Would it ever be useful to defer labeling?
           T1 #= T0 + T,
           Truth = true
         )
       ).

length_plus_1(Ls, MTD) :-
    length(Ls, MTDminus1),
    MTD #= MTDminus1 + 1.

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

state0_enrollment(Ls ^ Rs ^ Es, Ntotal) :-
    foldl(append, [Ls, Rs, Es], [], Cohorts),
    %% TODO: Maybe faster to sum each list, then sum the sums?
    maplist(\Q^N^(Q=_/N), Cohorts, Ns),
    sum(Ns, #=, Ntotal).

%?- ccd:state0_enrollment([] ^ [] ^ [], Ntot).
%@    Ntot = 0.
%?- ccd:state0_enrollment([1/6] ^ [1/4] ^ [2/5], Ntot).
%@    Ntot = 15.

ccd_state0_decision_state(CCD, Ls ^ [R | Rs] ^ Es, Decision, State) :-
    state0_enrollment(Ls ^ [R | Rs] ^ Es, Ntotal),
    enroll_max(Nmax),
    if_(Ntotal #< Nmax
	, ( tally_decision_ccd(R, Decision, CCD),
	    call(Decision, Ls ^ [R | Rs] ^ Es, State)
	  )
	, ( length_plus_1(Ls, MTD),
	    State = declare_mtd(MTD)
	  )
       ).

%% TODO: Consider restoring an 'mtd_notfound' concept to the trial.
%%       Even if this notion disappears 'WLOG' from a purely formal perspective,
%%       it remains part of the lingo, and with good reason. There really IS a
%%       difference between having found a moderate number of DLTs at the RP2D,
%%       and having probed nowhere into the toxic dose region.
%% TODO: But consider whether this type of linguistic expansion trespasses into a
%%       strictly *pharmacologic* realm which a formal analysis might best avoid.
%%       A good test may be to demonstrate that mtd_notfound truly adds something
%%       essential that can't be 'tacked on' in a post-processing step.
%%ccd_decisions(mtd_notfound(_)) --> [].
ccd_decisions(_, declare_mtd(_)) --> [].
ccd_decisions(CCD, S0) --> [(S0->A->S)],
			 { ccd_state0_decision_state(CCD, S0, A, S) },
			 ccd_decisions(CCD, S).

%% Examine the smallest possible trial -- a trial with just 1 dose!
%?- Trial+\(default_ccd(CCD), phrase(ccd_actions(CCD, []^[0/0]^[]), Trial)).
%@    Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[0/2]^[]),([]^[0/2]^[]->escalate->[]^[0/3]^[]),([]^[0/3]^[]->escalate->[]^[0/4]^[]),([]^[0/4]^[]->escalate->[]^[0/5]^[]),([]^[0/5]^[]->escalate->[]^[0/6]^[]),([]^[0/6]^[]->escalate->declare_mtd(1))]
%@ ;  Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[0/2]^[]),([]^[0/2]^[]->escalate->[]^[0/3]^[]),([]^[0/3]^[]->escalate->[]^[0/4]^[]),([]^[0/4]^[]->escalate->[]^[0/5]^[]),([]^[0/5]^[]->escalate->[]^[1/6]^[]),([]^[1/6]^[]->stay->declare_mtd(1))]
%@ ;  Trial = [([]^[0/0]^[]->stay->[]^[0/1]^[]),([]^[0/1]^[]->escalate->[]^[0/2]^[]),([]^[0/2]^[]->escalate->[]^[0/3]^[]),([]^[0/3]^[]->escalate->[]^[0/4]^[]),([]^[0/4]^[]->escalate->[]^[1/5]^[]),([]^[1/5]^[]->stay->[]^[1/6]^[]),([]^[1/6]^[]->stay->declare_mtd(1))]
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
state of ccd_decisions//1.

---

TODO: Note that 'cumulative-cohort matrices' (or 'cc matrices' for short) may
      better express the fact of our having eschewed a full path representation.
      OTOH, these matrices remain _summaries_ of the paths, and I'm not so sure
      that a term like 'path matrix' entails any kind of unique representation.

In terms of the actual tabular output to be read by R, a single row vector per
path seems ideal. The first element can be the selected dose, and this can be
followed by (T, N) pairs that even the human eye will be able to index easily.

*/

%% Rows of output table:
recdose_ccs, [] --> [(_->_->_^_^_)], recdose_ccs. % skip to end ...
recdose_ccs, [Pathvector] --> [ (S->_->declare_mtd(MTD)) ],
			      { recdose_state_pathvector(MTD, S, Pathvector) }.

recdose_state_pathvector(MTD, Ls^Rs^Es, Pathvector) :-
    reverse(Ls, Ks), % recall that Ls is a _descending_ list
    foldl(append, [Rs,Ks], Es, Qs),
    phrase(tallyvector(Qs), TNs),
    Pathvector = [MTD|TNs].

tallyvector([T/N|Qs]) --> [T, N], tallyvector(Qs).
tallyvector([]) --> [].

%% This predicate implements the common special case
%% where a CCD trial starts from the lowest dose.
ccd_d_pathvector(CCD, D, Pathvector) :-
    length(Tallies, D), maplist(=(0/0), Tallies),
    phrase(ccd_decisions(CCD, []^Tallies^[]), Path),
    phrase(recdose_ccs, Path, [Pathvector]).

ccd_d_cmax_nmax_tabfile(CCD, D, CohortMax, EnrollMax, Filename) :-
    atom_chars(File, Filename),
    format("Opening file ~q...~n", [File]), % feedback to console
    open(File, write, OS),
    format("Writing path vectors ..", []),
    (	retract(cohort_max(_)),
	fail
    ;	asserta(cohort_max(CohortMax))
    ),
    (	retract(enroll_max(_)),
	fail
    ;	asserta(enroll_max(EnrollMax))
    ),
    Ncols #= 2*D + 1,
    phrase(columns_format(Ncols), Format) ->
	'$cpu_now'(T0),
	(   ccd_d_pathvector(CCD, D, Pathvector),
	    format(OS, Format, Pathvector),
	    fail % exhaust all Pathvector solutions
	;   close(OS),
	    minutes_since(Minutes, T0),
	    format(".. done (~2f minutes).~n", [Minutes])
	).

% Copied from 'esc.pl'
columns_format(1) --> "~w~n".
columns_format(N) --> "~w\t",
		      { N #> 1,
			N1 #= N - 1 },
		      columns_format(N1).

:- op(900, xfx, ~>).

%% First, extract the path matrix:
path_matrix, [S ~> MTD] --> [ (S->_->declare_mtd(MTD)) ].
path_matrix, [] --> [(_->_->_^_^_)], path_matrix. % 'skip to end'

default_ccd(ccd([3/5, 4/8, 5/10, 6/12],		
		[1/1, 2/4, 3/6, 4/8, 5/11],
		[0/1, 1/8],
		12)).

%% I've implemented this predicate to more clearly exhibit the
%% sharing of arguments with ccd_state0_decision_state/4.
ccd_state0_matrix(CCD, State0, Matrix) :-
    phrase(ccd_decisions(CCD, State0), Path),
    phrase(path_matrix, Path, Matrix).

%% This predicate implements the common special case
%% where a CCD trial starts from the lowest dose.
ccd_d_matrix(CCD, D, Matrix) :-
    length(Tallies, D), maplist(=(0/0), Tallies),
    ccd_state0_matrix(CCD, []^Tallies^[], [Matrix]).

%?- Matrix+\(ccd:default_ccd(CCD), ccd_d_matrix(CCD, 2, Matrix)).
%@    Matrix = ([0/1]^[0/6]^[]~>2)
%@ ;  Matrix = ([0/1]^[1/6]^[]~>2)
%@ ;  Matrix = ([0/1]^[1/6]^[]~>2)
%@ ;  Matrix = ([0/1]^[2/6]^[]~>2)
%@ ;  Matrix = ([0/1]^[1/6]^[]~>2)
%@ ;  Matrix = ([0/1]^[2/6]^[]~>2)
%@ ;  Matrix = ([0/1]^[2/6]^[]~>2)
%@ ;  Matrix = ([]^[0/2,3/6]^[]~>1)
%@ ;  Matrix = ([]^[1/6,3/6]^[]~>1)
%@ ;  Matrix = ([]^[2/6,3/6]^[]~>1)
%@ ;  ...

%?- J+\(ccd:default_ccd(CCD), D=1, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 0.797 seconds
%@    % CPU time: 0.802 seconds
%@    J = 20.
%@    % CPU time: 0.549 seconds
%@    % CPU time: 0.554 seconds
%@    J = 20. % ^ PURE BASELINE
%@    % CPU time: 0.187 seconds
%@    % CPU time: 0.191 seconds
%@    J = 20.

%?- J+\(ccd:default_ccd(CCD), D=2, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 10.606 seconds
%@    % CPU time: 10.612 seconds
%@    J = 212.
%@    % CPU time: 5.878 seconds
%@    % CPU time: 5.882 seconds
%@    J = 212. % ^ PURE BASELINE
%@    % CPU time: 1.972 seconds
%@    % CPU time: 1.976 seconds
%@    J = 212.

%?- J+\(ccd:default_ccd(CCD), D=3, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 59.797 seconds
%@    % CPU time: 59.801 seconds
%@    J = 1151.
%@    % CPU time: 32.408 seconds
%@    % CPU time: 32.413 seconds
%@    J = 1151. % ^ PURE BASELINE
%@    % CPU time: 10.686 seconds
%@    % CPU time: 10.690 seconds
%@    J = 1151.

%?- J+\(ccd:default_ccd(CCD), D=4, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 388.514 seconds
%@    % CPU time: 388.518 seconds
%@    J = 6718.
%@    % CPU time: 62.763 seconds
%@    % CPU time: 62.767 seconds
%@    J = 6718.

%?- J+\(default_ccd(CCD), D=5, time(findall(M, ccd_d_matrix(CCD, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 370.651 seconds
%@    % CPU time: 370.655 seconds
%@    J = 39289.

regression :-
    default_ccd(CCD),
    J0s = [0, 20, 212, 1151, 6718, 26131], % 0-based list of values up to D=5
    D in 1..2, % 1..5,
    indomain(D),
    format(" D = ~d ...", [D]),
    time(findall(M, ccd_d_matrix(CCD, D, M), Ms)),
    length(Ms, J),
    format(" J(~d) = ~d.~n", [D,J]),
    nth0(D, J0s, J0),
    J #\= J0.

%?- asserta(ccd:cohort_max(6)), asserta(ccd:enroll_max(24)), ccd:regression.
%@  D = 1 ...   % CPU time: 1.089 seconds
%@    % CPU time: 1.093 seconds
%@  J(1) = 20.
%@  D = 2 ...   % CPU time: 13.192 seconds
%@    % CPU time: 13.196 seconds
%@  J(2) = 212.
%@ false.
