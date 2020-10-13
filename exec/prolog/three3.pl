% Attempt an enumeration of ALL 3+3 trials
:- use_module(library(clpfd)).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

/*

Aim initially to describe the 'Design 1' of Simon &al (1997).
The main difficulty in describing this design seems to be the
arbitrarily distant lookback to earlier cohorts, as required
to ensure that at most 1/6 toxicities occur at the declared
MTD.

What may suffice, however, is to define simple (Markov?) rules,
then PROVE that desired properties (such as '< 1/6') result.

ALTERNATIVELY (and surely *better*!), I might try to state the
desiderata, and allow the dose-escalation process to EMERGE as a
*result* of those specifications. For example, Skolnik &al (2008)
offer this definition:

> The MTD is defined as the dose level at which none or one
> of six participants (0% to 17%) experience a DLT, when at least
> two of three to six participants (33% to 67%) experience a DLT
> at the next highest dose.

To support such a declarative specification, I would have to
implement a scan of the full list generated thus far, perhaps
using semicontext notation as in the tree leaf-counting example
by Triska. The goals at any point would then be:
1. Can I declare an MTD?
2. If not, what should the next cohort be?
3. If I can't create another cohort, then declare "no MTD found".

NOTE HOW WONDERFULLY GOAL-DIRECTED THIS WOULD BE!!

Hmmm!! I think I begin to see implicit concepts REVEALED as I try
to think about the 3+3 in this manner:
a) A dose considered 'too toxic', initially this is (Hi+1)
b) A dose known to be 'safe', initially this is 0

Might I even get away with calling the 'safe' dose MTD?

REFERENCES

Skolnik, Jeffrey M., Jeffrey S. Barrett, Bhuvana Jayaraman, Dimple Patel,
and Peter C. Adamson. “Shortening the Timeline of Pediatric Phase I Trials:
The Rolling Six Design.” JCO 26, no. 2 (Jan 10, 2008): 190–95.
https://doi.org/10.1200/JCO.2007.12.7712.

=====

In all honesty, though, perhaps I should acknowledge that a
'running range' is carried forward implicitly as a state variable,
providing context for (de-)escalation decisions.

esc(D, Lo..Hi) --> ...

*/

% Let esc(D, Lo..Hi) *describe* a list of 3+3 cohorts starting from dose D in Lo..Hi.
% Read esc(D, Lo..Hi) as escalating FROM D to min(D+1,Hi).

%% TODO: WLOG set Lo == 1 by convention.

esc(Hi, Lo..Hi) --> [Hi * T], { Lo #< Hi, T in 0..3, indomain(T) }, % * is a ^ that went *splat*
		    (   { T in 0..1 } -> [mtd_notfound(Hi)]
		    ;   des(Hi, Lo)
		    ).
esc(D, Lo..Hi) --> [D1 ^ T], { D1 #= D + 1, D1 in Lo..Hi, T in 0..3, indomain(T) },
		   (   { T #= 0 } -> esc(D1, Lo..Hi)
		   ;   { T #= 1 } -> sta(D1, Lo..Hi)
		   ;   des(D1, Lo)
		   ).

% TODO: Do I really need the special case of sta(D, _..D)?
sta(D, _..D) --> [D - 0], [mtd_notfound(D)].
sta(D, Lo..Hi) --> [D - 0], { D #< Hi, D in Lo..Hi, indomain(D) }, esc(D, D..Hi).
sta(D, Lo.._) --> [D - T], { T in 1..3, indomain(T) }, des(D, Lo).

% As a mirror image of esc, des(D, Lo) moves downward FROM D to max(D-1,Lo).
% NB: Implicit in de-escalation to D-1 is that Hi #= D - 1.
des(D, Lo) --> { D_1 #= D - 1 },
	       (   { D_1 #= Lo } -> [declare_mtd(Lo)]
	       ;   [D_1 : T], { T in 0..3, indomain(T) },
		   (   { T in 0..1 } -> [declare_mtd(D_1)]
		   ;   des(D_1, Lo)
		   )
	       ).

%% TODO: Write cohorts as (Dose ^ Tox / Enrolled).
%%       This will enable cohorts of size other than 3,
%%       as needed e.g. to describe an accelerated titration phase.

% n_trials(+Drange, XY)
n_trials(Drange, XY) :-
    Dmax in Drange, indomain(Dmax),
    findall(Tr, phrase(esc(0, 0..Dmax), Tr), Trials),
    length(Trials, N),
    XY = (Dmax, N).

/*******

I don't feel entirely happy with the way I have intertwined the
escalation rules with MTD-declaration logic above. It would be
preferable (if feasible) to define alternating 'coroutines' of
dose escalation and attempted MTD declaration.

Under such a scheme, we would at each step attempt to declare an
MTD, and -- failing that -- apply the dose-escalation rules.

A quite natural way to implement this, perhaps, is to build the
list of cohorts 'in reverse', with the latest cohort at the head
of the list.

IDEALLY, the program would directly express the *intuition* beneath
the dose-escalation rules, rather than merely implement the rules.

Another interesting approach would be to try to PROVE various local
properties of the existing DCG.

*/

de([], [1^T]) :- T in 0..3, indomain(T).
de([D^0|E], [D1^T,D^0|E]) :- D1 #= D + 1, D1 #=< 4, T in 0..3, indomain(T).
de([D^1|E], [D^T,D^1|E]) :- T in 0..3, indomain(T).
de([D^T|E], [D_1^T_1,D^T|E]) :- T in 2..3, D_1 #= D - 1, D_1 #> 0, T_1 in 0..3, indomain(T_1).

/*******

What if I acknowledge the dose-centeredness of these designs,
and build my program on this concept?

The representation of dose range becomes important.

One simple way to do it would be to write a list initially of the form
[0/0, 0/0, ..., 0/0], with one slot per nonzero dose level.

We then reason over this representation, implementing state transformations
according to the pharmacologic intuitions behind 'judging' dose levels by
the counts of toxicities.

What categorical 'judgements' do we make?

a. A dose level with 0/6 or 1/6 toxicities is SAFE
b. A dose level with 0/3 toxicities is LOW -- i.e., too low to enroll next cohort
c. A dose level with 2+ toxicities is UNSAFE

I believe that these are the underlying pharmacologic intuitions of the 3+3 &al.

The basic attitude--or *approach*--of 3+3 is to whittle away at our uncertainty
about some prespecified set of dose levels, until it is fully resolved.

We recognize that resolution in the form of:
1. A list [SAFE, UNSAFE] in which case the SAFE dose in MTD ** def'n of MTD! **
2. A singleton list with SAFE dose => no MTD found
3. A singleton list with UNSAFE dose => MTD is ZERO

Is there any advantage to seeing the list as a STACK? One problem with this is
that it's likely to impinge on generality. Nevertheless, it does nicely express
the intuition that we work stepwise from the 'safe end' of the dose range, which
is foundational to the dose-escalation concept. Thus, 'entangling' the program
implementation with this intuition could even be a good thing if 'speculative
generality' continues to apply as a valid criticism in logic programming!

*/

%% Hmmm... Maybe I need to reference higher doses here, as well!
%% Perhaps a dose is safe also if any HIGHER dose is safe.
:- true.

% Notice how nicely the code below lays out our dose-finding intuitions!
% This is EXACTLY the type of benefit one hopes to gain from Prolog.
% Observe how all the reasoning happens over ASCENDING SEQUENCES OF DOSES;
% this is an essential characteristic of the underlying intuitions.

% Note that we would avoid the need for some special-case handling
% if we imposed the interpretation on the 'missing' dose past the end
% of the list, that it represents a dose known to be unsafe.
% But this need not hold GENERALLY, since another reason for capping
% the dose sequence considered is that any higher doses are known to
% be on an efficacy plateau. Thus, this *special* interpretation, tho
% leading to 'more elegant' code, harms the generality of the program.

safe(T/6) :- T #=< 1. % TODO: Is this 'defaulty'?
safe([D|H]) :-  % Head dose in a sequence is safe, provided
    (	safe(D) % either D < 2/6 ..
    ;	safe(H) % ..or some later dose meets that criterion.
    ).
unsafe(T / N) :- N #>= 3, T #> 1. % 2+/_ is unsafe
% The head dose level in a sequence is too low to enroll if
low([D|_]) :- safe(D). % ..it's safe IN ITSELF,
low([0/3, T/N | _]) :- \+ unsafe(T/N). % ..or 0/3 and next-higher dose is not unsafe.
% Note in particular that 0/3 without a successor (i.e., 0/3 at top dose)
% is NOT too low to enroll. This little predicate does a lot of work,
% since it's what forces escalation.

mtd([S,U|_]) :- % MTD is found at head of dose sequence considered
    safe(S),
    unsafe(U).
mtd([U|_], 0) :- unsafe(U).
mtd([S,U|_], 1) :- mtd([S,U|_]).
mtd([_|E], N) :-
    N #> 1,
    N_1 #= N - 1,
    mtd(E, N_1).

% Enrolling N and observing T toxicities updates our dose range
% from R0 to R, and said enrollment was indeed permitted by R0.
%%range0_cohort_range(R0, T/+N, R)

range0_cohort_range([L|R0], T/N, [L|R]) :- % we bypass low dose levels L
    low([L|R0]),
    range0_cohort_range(R0, T/N, R).

range0_cohort_range([T0/N0 | R0], T/N, [T1/N1 | R0]) :-
    \+ low([T0/N0 | R0]), % why doesn't failure of above clause ensure this?
    \+ unsafe(T0/N0),
    T in 0..N, indomain(T),
    T1 #= T0 + T,
    N1 #= N0 + N.

% Now define a trial RECURSIVELY, with an eye to employing DCG
% for state threading ...

% A completed trial is a list of dose ranges with successive updates
% consistent with the dose-escalation rules.
initialize(1, [0/0]).
initialize(N, [0/0 | S]) :- N #> 1, N_1 #= N - 1, initialize(N_1, S).
trial(N) --> { initialize(N, Start) }, [Start], step(Start).
step(R0) -->
    (	{ mtd(R0, MTD) } -> [declare_mtd(MTD)]
    ;	{ range0_cohort_range(R0, _/3, R) }, [R], step(R) % preserve choice point!
    ;	{ \+ range0_cohort_range(R0, _/3, _) } -> [mtd_notfound(D)], { length(R0, D) }
    ).

% I do like the pattern of describing a step via range0_cohort_range,
% and then leaving mtd_notfound as a kind of 'failure mode'. But the
% SYNTAX of the above feels deficient, or at least 'inelegant'!

n_trials_both(Drange, XAB) :-
    Dmax in Drange, indomain(Dmax),
    findall(Tr, phrase(esc(0, 0..Dmax), Tr), TrialsA),
    length(TrialsA, Na),
    findall(Tr, phrase(trial(Dmax), Tr), TrialsB),
    length(TrialsB, Nb),
    XAB = (Dmax, Na, Nb).


%% Transform the A path representations to the arrays T(c,d,j).

path_matrix(P, D, M) :-
    length(C1,D),
    length(C2,D),
    M = (C1,C2),
    path_matrix_(P, M),
    maplist(ground_or_nil, C1),
    maplist(ground_or_nil, C2).

path_matrix_([D^T|Rest], (C1,C2)) :-
    nth1(D, C1, T),
    path_matrix_(Rest, (C1,C2)).

path_matrix_([D-T|Rest], (C1,C2)) :-
    nth1(D, C2, T),
    path_matrix_(Rest, (C1,C2)).

path_matrix_([D*T|Rest], (C1,C2)) :- path_matrix_([D-T|Rest], (C1,C2)).
path_matrix_([D:T|Rest], (C1,C2)) :- path_matrix_([D-T|Rest], (C1,C2)).

path_matrix_([], (_,_)).

ground_or_nil(Term) :- ground(Term).
ground_or_nil(nil) :- true.

%% Write out separate R input files for D in 2..7
