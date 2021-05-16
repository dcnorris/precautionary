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
:- use_module(library(pairs)).
:- use_module(library(dcgs)).
:- use_module(library(lambda)).
:- use_module(library(time)).
:- use_module(library(debug)).
:- use_module(library(between)).
:- use_module(library(dif)).
:- use_module(library(reif)).

% TODO: Review library(debug); it includes * and other useful predicates

%% -----------------------------------------------------------------------------

% My initial emphasis is on generating all possible paths (CPE) for the BOIN
% design set forth in the table above. Although the BOIN design of [1] lacks
% any terminating principle except elimination of all doses, we do need such
% a rule here. The most natural type of rule, in view of the Table above,
% might be a 'stop-for-consensus' type of rule as found in package 'dtpcrm'.
% This is specified as a maximum number of patients to enroll at any 1 dose.

%% In general, I may require a max-enrollment *parameter*,
%% however unsightly it may be tagging along like this ...
tally_maxn(DLTs/Enrolled, MaxN) :-
    ground(MaxN),
    Enrolled in 0..MaxN,
    indomain(Enrolled),
    DLTs in 0..Enrolled,
    indomain(DLTs).

%% Having a 'default' max cohort size will make dev & test a bit easier:
tally(DLTs/Enrolled) :- tally_maxn(DLTs/Enrolled, 12).

%% tally/1 terminates:
%?- tally(_), false.
%@ false.

%% -----------------------------------------------------------------------------

/*

I must define the comparison between a tally and a boundary, so that all
possible true cases are taken into account. Mathematically, this means
a complete mapping from {tallies} x {boundaries} --> {<,=,>}, much like
the one achieved by zcompare/3.

Logically, defining the set B = {boundaries} is a prerequisite to this!

Perhaps I even need to examine B as a *group* with binary operations
including min() and max(). Suitable definitions here might enable me to
describe 'minimal' representations and equivalence of b1, b2 ϵ B.

FIRSTLY, the 'or' compositionality of the boundary lists needs to be
explored. Perhaps these boundary comparisons are just natural extensions
of the earlier comparisons between individual tallies. (OTOH, maybe this
extension is not simple, if it yields a more complete comparison relation.)

SECONDLY, I should distinguish between the full representation of a tox
boundary, and various *abbreviations*. Perhaps the boundaries *formally*
are infinite lists of numerators over ℕ, and perhaps any such abbreviation
ought to be (until we worry about performance) expanded to this CANONICAL
form before any computations are done.

*/

%% 1. TODO: Define the CANONICAL FORM of tox boundaries
/*

Ideally, a tox boundary will be comparable with ANY given tally.
Thus, tox boundaries eliminate the incomparabilities that in general
confound comparisons (on the reachability principle) between individual
tallies having different denominators.

This is done by eliminating 'gaps in the fenceposts' of the boundaries.

OOH! Maybe a boundary is (at least in the *abstract*) best understood
as a function that yields a tally numerator for any denominator.
That is, the boundary provides a 'common-denominator' operation.

Not every list of the form [_/0, _/1, _/2, _/3, ...] will necessarily
be a *useful* one, but perhaps they are all 'valid' boundaries, under
the most elegant description. The list just becomes a compendium of all
possible comparisons, generating a complete function {tallies} --> {<,=,>}.

A crucial question of SEMANTICS (i.e., the *meaning* of 'boundary') can
be dealt with by asking which such lists are sensible/useful.

It doesn't hurt to _begin_ thinking about the abbreviation, even at this
early stage. What would the singleton list [0/2] generate by way of an
infinite list of the form [_/0, _/1, _/2, _/3, ...]? Might it generate
different lists, depending on the comparison (</=/>) being attempted?

Let me think SPATIALLY about this! How does a single tally, plotted in
the (T,N) plane, divide that plane up?

*/

qcompare(~~, T1/N1, T2/N2) :- % REACHABILITY
    tally(T1/N1), tally(T2/N2),
    (	N1 #= N2,
	T1 #= T2
    ;	N1 #< N2,
	T1 #=< T2, T2 #=< T1 + (N2 - N1)
    ;	N1 #> N2,
	T1 #>= T2, T2 #>= T1 + (N2 - N1)
    ).

qcompare(=<, T1/N1, T2/N2) :-
    tally(T1/N1), tally(T2/N2),
    (	N1 #>= N2, T1 #=< T2
    ;	N1 #< N2, N1 - T1 #>= N2 - T2
    ).
	
qcompare(<, T1/N1, T2/N2) :-
    tally(T1/N1), tally(T2/N2),
    (	N1 #>= N2, T1 #< T2
    ;	N1 #< N2, N1 - T1 #> N2 - T2
    ).
	
qcompare(>=, T1/N1, T2/N2) :-
    tally(T1/N1), tally(T2/N2),
    (	N1 #=< N2, T1 #>= T2
    ;	N1 #> N2, N1 - T1 #=< N2 - T2
    ).
	
qcompare(>, T1/N1, T2/N2) :-
    tally(T1/N1), tally(T2/N2),
    (	N1 #=< N2, T1 #> T2
    ;	N1 #> N2, N1 - T1 #< N2 - T2
    ).

%% Demonstrate that strict inequalities are exclusive of ~~
%?- time((qcompare(~~, Q1, Q2), qcompare(<, Q1, Q2))).
%@    % CPU time: 37.759 seconds
%@ false.
%?- time((qcompare(~~, Q1, Q2), qcompare(>, Q1, Q2))).
%@    % CPU time: 37.625 seconds
%@ false.

%% Show that =< and >= hold simultaneously only in case of equivalence:
%?- time((tally(Q2), qcompare(>=, Q1, Q2), qcompare(=<, Q1, Q2), dif(Q1, Q2))).
%@    % CPU time: 56.061 seconds
%@ false.

/*

WHERE SENSIBLE, extend the above comparisons to LISTS of arg 2 ...

What is the right semantics for qcompare(C, Q, [Q1, Q2])?
Obvious candidates are (qcompare(C, Q, Q1), qcompare(C, Q, Q2)) ("and")
and (qcompare(C, Q, Q1); qcompare(C, Q, Q2)) -- ("or").

Certainly, the "or" case is more 'efficient' (because *lazy*).

But what do the BOIN tox boundaries require, 'as written'?

Ahh, suppose I wanted to allow the tox boundary to be 'staked out'
with one tally T/N per N ∊ ℕ. In that case, it would suffice to demonstrate
the inequality for just the matching element. This naturally points toward
the "or" connectivity. In this case, the implementation is straightforward
as well, in terms of lists.

But might the strict and non-strict inequalities naturally require opposite
semantics? If 'hitting a boundary' is something that happens if it happens
at any point along the boundary, then NOT hitting the boundary seems to
require a logical CONJUNCTION.

*/

%% Extend this relation to a list arg #3. Note that the strict inequalities
%% demand a logical CONJUNCTION, consistent with 'not hitting the boundary'
%% demanding that we miss EVERY point 'staked out' along the boundary.
%% Conversely, 'hitting the boundary' means hitting ANY point.
qcompare(<, Q1, []) :- tally(Q1).
qcompare(<, Q1, [Q | Qs]) :-
    qcompare(<, Q1, Q),
    qcompare(Q1, Qs).

qcompare(>, Q1, []) :- tally(Q1).
qcompare(>, Q1, [Q | Qs]) :-
    qcompare(>, Q1, Q),
    qcompare(Q1, Qs).

%?- if_(=(1, 1), A = 1, A = 0).
%@    A = 1.
%?- if_(1 = 0, A=1, A=0).
%@    A = 0. % Dang!

%% TODO: Here is where (apparently) I would use if_/3.
%%       It seems this would require reifying the truth-value
%%       of qcompare/3 in a new, 4th argument. But why wouldn't
%%       this amount to the same work as =/3 in reif.pl, which
%%       itself uses (->) twice?
qcompare(=<, Q1, [Q]) :- qcompare(=<, Q1, Q).
qcompare(=<, Q1, [Q | Qs]) :-
    (	qcompare(=<, Q1, Q) -> true
    ;	qcompare(=<, Q1, Qs)
    ).

qcompare(>=, Q1, [Q]) :- qcompare(>=, Q1, Q).
qcompare(>=, Q1, [Q | Qs]) :-
    (	qcompare(>=, Q1, Q) -> true
    ;	qcompare(>=, Q1, Qs)
    ).

/*

Interestingly, the analysis via qcompare/3 (and the proofs it enables)
assist greatly in the *analysi* of the problem, but probably do not
enter into the final *representation*.

Rather, the representation ought to exploit the understanding emerging
from the analysis that led up to qcompare/3, and exploit it in such a
way that it 'factors out', becoming invisible.

This 'disappearing act' happens because the new understanding supports
the implicit (thus, invisible) claim that the representation is sufficient.

I will claim that all reasonable decision boundaries (whether floors or
ceilings) look like ascending stairs with step heights of 1. This permits
me to represent any [reasonable] boundary by a list giving the positions
of the steps. The semantics of floors and ceilings are even mirror images:

A FLOOR represented as [1, 8] expands to 'canonical form' as follows:

[1, 8] => [0/1, 1/8] => [0/1, ..., 0/7, 1/8, 1/9, ..., 1/Cmax].

Using the analysis embodied in qcompare/3, we could demonstrate that

Q =< [0/1, 1/8] iff Q =< [0/1,...,0/7,1/8,...].

This is simply by transitivity:

Q =< 0/1 ==> Q =< 0/N =< 0/1 for N > 1;
Q =< 1/8 ==> Q =< 1/N =< 1/8 for N > 8.

A CEILING represented as [3, 6, 9, 11] expands to canonical form like so:

[3, 6, 9, 11] => [3/3, 4/6, 5/9, 6/11]
              => [3/3, 3/4, 3/5,
                  4/6, 4/7, 4/8,
                  5/9, 5/10,
                  6/11, ..., 6/Cmax].

The more direct expression of this logic is that the first entry T/N
in a ceiling list is WLOG on the worst-case diagonal T==N. Thus, the
first entry in the list is both T and N to start. As with a tox floor,
the T's ascend in steps of 1.

T = [3, 4, 5, 6]
N = [3, 6, 9, 11]

=> [3/3, 4/6, 5/9, 6/11] => ...

If it were desired to obtain a comparator list suitable for use as arg 3
in qcompare/3, then this could be done by a further step:

=> [3/5, 4/8, 5/10, 6/Cmax].

But this is probably not helpful. Under a "state-what-holds" approach
to the programming, we want instead simple CLP(ℤ) constraints to define
the dose-escalation decisions.

... UNLESS(!!) if_/3 comes to the rescue...

(If this turns out to work well, then in this very practical application
of Prolog, we have a case where using the if_/3 predicate preserves a
certain type of 'thoroughbass' of formality. TODO: Revisit this intuition
at some later point.)

*/

hit_ceiling_t(Q, [], false).
hit_ceiling_t(Q, [C|Cs], Truth) :-
    (	qcompare(>=, Q, C) -> Truth = true
    ;	hit_ceiling_t(Q, Cs, Truth)
    ).

hit_floor_t(Q, [], false).
hit_floor_t(Q, [F|Fs], Truth) :-
    (	qcompare(=<, Q, F) -> Truth = true
    ;	hit_floor_t(Q, Fs, Truth)
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

%% Find concise expressions of floor and ceiling tally lists
/*
In a floor-type boundary list, if there are 2 elements Q1 =< Q2,
then Q1 is superfluous, since it is 'under the floor', and the
transitivity of (=<) [this follows from the CONVEXITY of the
sub-lattices it generates] means that Q2 'has Q1 covered'.
*/
concise_floor(Concise, Floor) :-
    zip(_, Ns, Floor), % Ns are denominators of Floor
    tallies_sorted(Floor, FloorSorted),
    concise_floor_(Concise, FloorSorted).

%% TODO: Actually implement the reduction!
%%       What's needed is some kind of reduction at the head,
%%       which because of the sort ought to be possible in a
%%       single linear pass.
concise_floor_([_], [_]).
%concise_floor_([Q1|Qs], [Q1,Q2|Qs]) :- ... % TODO

tallies_sorted(Qs, Sorted) :-
    zip(_, Ns, Qs),
    %% TODO: Could use map_list_to_pairs/3 here instead,
    %%       provided I defined a den/2 predicate or that
    %%       2nd-arg-of-functor / is already available.
    pairs_keys_values(NQs, Ns, Qs),
    sort(NQs, NSorted),
    pairs_keys_values(NSorted, _, Sorted).

%?- tallies_sorted([0/7, 0/4, 0/5, 0/1, 0/3, 0/6, 0/2, 1/10, 1/9, 1/8, 1/12, 1/11], S).
%@    S = [0/1,0/2,0/3,0/4,0/5,0/6,0/7,1/8,1/9,... / ...|...]
%@ ;  false.

%% 2. TODO: Define ABBREVIATIONS of tox boundaries via abbrev_boundary/2

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
%   T/N &=< [1, 8, 15],
% although I am inclined to reject the latter as more 'cute' than *clear*.
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
    chain(#<, Ns),
    maplist(#>(0), Ns),
    %% Note there are further constraints on what might be considered
    %% reasonable tox boundaries. There shouldn't be duplicated Ns, e.g.,
    %% and all the Ts should be non-negative. But perhaps these are best
    %% postponed until some kind of *search* over the space of possible
    %% tox boundaries becomes desirable.
    maplist(#=<, Ts, Ns).

%?- sort([A,B,C], [3,1,2]).
%@    A = 3, B = 1, C = 2.

%?- sort([A,B,C], [3,1,2]), sort([A,B,C], [3,1,2]).
%@ false.

%?- tox_boundary([0/2, 1/4]).
%@ false.
%@    true
%@ ;  false.

%?- tox_boundary([0/5, 1/4]).
%@ false.
%@ false.

%% TODO: Why does sortedness of Ns not appear in constraints below?
%%       (Do I need to assert sortedness in terms clpz that understands?)
%?- tox_boundary(Bdy).
%@    Bdy = []
%@ ;  Bdy = [_B/_A], clpz:(_A#>=_B), clpz:(_A in inf.. -1), clpz:(_B in inf.. -1)
%@ ;  Bdy = [_C/_A,_D/_B], clpz:(_A#=<_B+ -1), clpz:(_A#>=_C), clpz:(_B#>=_D), clpz:(_A in inf.. -2), clpz:(_C in inf.. -2), clpz:(_B in inf.. -1), clpz:(_D in inf.. -1)
%@ ;  ...
%@    Bdy = []
%@ ;  Bdy = [_B/_A], clpz:(_A#>=_B), clpz:(_A in inf.. -1), clpz:(_B in inf.. -1)
%@ ;  Bdy = [_B/_A,_D/_C], clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_A in inf.. -1), clpz:(_B in inf.. -1), clpz:(_C in inf.. -1), clpz:(_D in inf.. -1)
%@ ;  Bdy = [_B/_A,_D/_C,_F/_E], clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_E#>=_F), clpz:(_A in inf.. -1), clpz:(_B in inf.. -1), clpz:(_C in inf.. -1), clpz:(_D in inf.. -1), clpz:(_E in inf.. -1), clpz:(_F in inf.. -1)
%@ ;  ...
%@    Bdy = []
%@ ;  Bdy = [_B/_A], clpz:(_A#>=_B), clpz:(_A in inf.. -1), clpz:(_B in inf.. -1)
%@ ;  Bdy = [_B/_A,_D/_C], clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_A in inf.. -1), clpz:(_B in inf.. -1), clpz:(_C in inf.. -1), clpz:(_D in inf.. -1)
%@ ;  ...
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

escalate(Q) :- hit_floor(Q, [0/1, 1/8]).
deescalate(Q) :- hit_ceiling(Q, [1/1, 2/4, 3/6, 4/8, 5/10, 6/12]).
remove(Q) :- hit_ceiling(Q, [3/5, 4/8, 5/10, 6/12]).

%% TODO: Note the duplication between remove/2 and deescalate/2.
%%       Consider whether these should be refined to avoid this.
%% TODO: Might a zcompare/3-like idiom further the cause of purity?

/*
Even as I set out to write the first lines of this predicate,
I can see the value of having complementary &** predicates
that enable me to state definitively whether a boundary has
been breached, or not.
In a way, what I need is a kind of zcompare/3 for tox boundaries!
Let me not impose burdens on state-machine code, that more properly
belong to basic comparison operators!
I don't think I have yet fully worked out the MEANING of toxicity
boundaries as used in BOIN.
Keep in mind that I want the code ultimately to be verifiable
on inspection by Ying Yuan and others -- at least where it is not
self-verifying by proofs executed in Prolog itself.
*/

%% NB: The left-hand list of Ls ^ Rs is sorted in descending order.
%%     So heads L and R in [L|Ls] ^ [R|Rs] are tallies belonging to
%%     *adjacent* doses, notwithstanding their non-juxtaposition in
%%     our left-right reading of the term.

cohort_full(N, true) :- N #>= 12.
cohort_full(N, false) :- N in 0..11.
cohort_full(_/N, TF) :- cohort_full(N, TF).

enroll(T0/N0, T1/N1) :-
    cohort_full(N0, false),
    N1 #= N0 + 1,
    T in 0..1, % T is the pending tox assessment of newly-enrolled patient
    T1 #= T0 + T.

%% Overload enroll/2 on lists of tallies, enrolling head tally
enroll([T0/N0 | Qs], [T1/N1 | Qs]) :-
    enroll(T0/N0, T1/N1).

length_plus_1(Ls, MTD) :-
    length(Ls, MTD_1),
    MTD #= MTD_1 + 1.

stay(Ls ^ [R | Rs], Ls ^ [R1 | Rs]) :- enroll(R, R1).
stay(Ls ^ [_/N | _], declare_mtd(MTD)) :- cohort_full(N, true),
					  length_plus_1(Ls, MTD).

escalate(Ls ^ [T/N], State) :- % NB: this is a 'clamped' situation
    stay(Ls ^ [T/N], State).

escalate(Ls ^ [Q, R | _], [Q | Ls] ^ [R1 | _]) :- enroll(R, R1).
escalate(Ls ^ [Q, R | _], declare_mtd(MTD)) :- cohort_full(R, true),
					       length_plus_1([Q|Ls], MTD).

deescalate([] ^ _, declare_mtd(0)). % deescalate from already-lowest dose

%% TODO: Use 'guard' to parse this into separate clauses.
deescalate([L | Ls] ^ Rs, Ls ^ [L1 | Rs]) :- enroll(L, L1).
deescalate([L | Ls] ^ _, declare_mtd(MTD)) :- cohort_full(L, true),
					      length(Ls, MTD).

%% TODO: Note how TRIVIAL state0_action_state/3 has become.
%%       Does this demand refactoring, or does it represent other
%%       types of opportunity -- even meta-interpretation?

state0_action_state(Ls ^ [R | Rs], escalate, State) :-
    escalate(R),
    escalate(Ls ^ [R | Rs], State).

state0_action_state(Ls ^ [R | Rs], deescalate, State) :-
    deescalate(R),
    deescalate(Ls ^ [R | Rs], State).

state0_action_state(Ls ^ [R | Rs], stay, Ls ^ [R1 | Rs]) :-
    %% This is not sound if R is not ground!
    %% zcompare/3 caters to all possible case, and never loses
    %% any of them or incorrectly commits to them.
    %% Always have all correct cases in mind.
    %% You can see by looking at this code that it is incomplete!
    %% Use clauses to enumerate the possible cases.
    %% Unsafe to use extralogical constructs (e.g., sort/2).
    %% TODO: Look at term rewriting video (where if-then-else used)
    %% Consider if_ here also.
    %% But no matter what, be able to state precisely what holds;
    %% a linguistic task: formulate precisely what conditions
    %% make this true.
    %% TODO: Use the MGQ to find where this unsound clause fails!
    (	escalate(R) -> false
    ;	deescalate(R) -> false
    ;	%% Until I gain greater (e.g., zcompare/3-style) control over the
	%% boundary-hitting concept, I must treat 'stay' as a DEFAULT action:
	enroll(R, R1)
    ).

%% Note how this state-machine naturally starts up from a blank slate:
%?- state0_action_state([] ^ [0/0, 0/0, 0/0], Action, State).
%@    Action = stay, State = []^[_A/1,0/0,0/0], clpz:(_A in 0..1)
%@ ;  false.

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

%% Examine the smallest possible trial -- a trial with just 1 dose!
%?- phrase(actions([]^[0/0]), Trial).
%@    Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  ...
%@    Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  ...
%@    Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
%@ ;  Trial = [([]^[0/0]->stay->[]^[0/1]),([]^[0/1]->escalate->[]^[0/2]),([]^[0/2]->escalate->[]^[0/3]),([]^[0/3]->escalate->[]^[0/4]),([]^[0/4]->escalate->[]^[0/5]),([]^[0/5]->escalate->[]^[0/6]),([]^[0/6]->escalate->[]^[0/7]),([]^[... / ...]->escalate->[]^[...]),([]^ ... ->escalate-> ... ^ ...),(... -> ...)|...]
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
