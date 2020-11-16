% Attempt an enumeration of ALL 3+3 trials
:- use_module(library(clpfd)).
:- use_module(library(pio)).

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

tox(T) :- T in 0..3,
	  *indomain(T). % * switches off labeling

% Mnemonic: * is ^ that 'splatted' on dose ceiling.
esc(Hi, Lo..Hi) --> [Hi * T], { Lo #< Hi,
				tox(T) },
		    (  {T #=< 1}, [mtd_notfound(Hi)]
		    ;  {T #>= 2}, des(Hi, Lo)
		    ).
esc(D, Lo..Hi) --> [D1 ^ T], { D1 #= D + 1,
			       D1 in Lo..Hi,
			       tox(T) },
		   (  {T #= 0}, esc(D1, Lo..Hi)
		   ;  {T #= 1}, sta(D1, Lo..Hi)
		   ;  {T #> 1}, des(D1, Lo)
		   ).

% TODO: Is special case sta(D, _..D) needed?
sta(D, _..D) --> [D - 0], [mtd_notfound(D)].
sta(D, Lo..Hi) --> [D - 0], { D #< Hi,
			      D in Lo..Hi },
		   esc(D, D..Hi).
sta(D, Lo.._) --> [D - T], { tox(T), T #> 0 },
		  des(D, Lo).

% As a mirror image of esc//2, des(D, Lo) moves
% downward FROM D, to max(D-1,Lo).
% NB: De-escalation to D-1 implies Hi #= D - 1.
%% TODO: Does this mean des//2 could be written as des(Lo..D),
%%       somehow exploiting the impossibility of Y..X with Y > X?
des(D, Lo) --> { D_1 #= D - 1 },
	       (  {D_1 #= Lo}, [declare_mtd(Lo)]
	       ;  {D_1 #> Lo}, [D_1 : T], {tox(T)},
		  (  {T #=< 1}, [declare_mtd(D_1)]
		  ;  {T #>= 2}, des(D_1, Lo)
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

?- n_trials(1..8, XY).
%@ XY =  (1, 5) ;
%@ XY =  (2, 15) ;
%@ XY =  (3, 37) ;
%@ XY =  (4, 83) ;
%@ XY =  (5, 177) ;
%@ XY =  (6, 367) ;
%@ XY =  (7, 749) ;
%@ XY =  (8, 1515).
%@ XY =  (1, 5) ;
%@ XY =  (2, 15) ;
%@ XY =  (3, 37) ;
%@ XY =  (4, 83) ;
%@ XY =  (5, 177) ;
%@ XY =  (6, 367) ;
%@ XY =  (7, 749) ;
%@ XY =  (8, 1515).
%@ XY =  (1, 10) ;
%@ XY =  (2, 46) ;
%@ XY =  (3, 154) ;
%@ XY =  (4, 442) ;
%@ XY =  (5, 1162) ;
%@ XY =  (6, 2890) ;
%@ XY =  (7, 6922) ;
%@ XY =  (8, 16138).

?- phrase(esc(0, 0..2), Tr).
%@ Tr = [1^0, 2^0, 2*_7658, mtd_notfound(2)],
%@ _7658 in 0..1 ;
%@ Tr = [1^0, 2^0, 2*_9396, 1:_9408, declare_mtd(1)],
%@ _9396 in 2..3,
%@ _9408 in 0..1 ;
%@ Tr = [1^0, 2^0, 2*_10960, 1:_10972, declare_mtd(0)],
%@ _10960 in 2..3,
%@ _10972 in 2..3 ;
%@ Tr = [1^0, 2^1, 2-0, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^1, 2-_13576, 1:_13588, declare_mtd(1)],
%@ _13576 in 1..3,
%@ _13588 in 0..1 ;
%@ Tr = [1^0, 2^1, 2-_15140, 1:_15152, declare_mtd(0)],
%@ _15140 in 1..3,
%@ _15152 in 2..3 ;
%@ Tr = [1^0, 2^_17066, 1:_17078, declare_mtd(1)],
%@ _17066 in 2..3,
%@ _17078 in 0..1 ;
%@ Tr = [1^0, 2^_18618, 1:_18630, declare_mtd(0)],
%@ _18618 in 2..3,
%@ _18630 in 2..3 ;
%@ Tr = [1^1, 1-0, 2^0, 2*_20676, mtd_notfound(2)],
%@ _20676 in 0..1 ;
%@ Tr = [1^1, 1-0, 2^0, 2*_21998, declare_mtd(1)],
%@ _21998 in 2..3 ;
%@ Tr = [1^1, 1-0, 2^1, 2-0, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-_23998, declare_mtd(1)],
%@ _23998 in 1..3 ;
%@ Tr = [1^1, 1-0, 2^_25308, declare_mtd(1)],
%@ _25308 in 2..3 ;
%@ Tr = [1^1, 1-_26568, declare_mtd(0)],
%@ _26568 in 1..3 ;
%@ Tr = [1^_27854, declare_mtd(0)],
%@ _27854 in 2..3.

%@ XY =  (1, 10) ;
%@ XY =  (2, 46) ;
%@ XY =  (3, 154) ;
%@ XY =  (4, 442) ;
%@ XY =  (5, 1162) ;
%@ XY =  (6, 2890) ;
%@ XY =  (7, 6922) ;
%@ XY =  (8, 16138). 

% The DCG esc//2 is a dose-escalation trial proceeding
% from a condition where the current dose in arg1, and
% the dose range still being entertained as a possible
% MTD is arg2.

% Had I tried to write this less 'economically', I may
% have preferred something like:
% below_current_above//3
% below_current_above(Ls, D, Rs)
%% We start a dose-escalation trial by escalating from the
%% null situation of 0 toxicities at dose zero ([]).
detrial(N) --> { length(Ds, N) }, de(>, Ds, [], []).
% de(C, D, E) means we just observed a tox count (<,=,>)~(low,middling,high)
% on dose length(D), with additional doses E above it still being considered.
% Perhaps the full computation requires maintaining 3 lists, the first of which
% is doses known to be safe, the second being those tried only once (and therefore
% not presumed safe or unsafe) and the third being doses not yet tried.
% That first list accomplishes essentially the tracking of the lower end
% of the currently considered range.
% Could I aim for greater clarity on this?
% What if I think about the middle list as a temporary holding area for doses
% that have been tried only once, with [in 3+3, at least] no toxicities?
% Can I treat them initially in some 3+3-specific way, but with the hope of
% finding more general treatments later?
% LIST 1: Descending list of doses with head having 0/6 or 1/6
% LIST 2: Descending list of doses all having resulted so far 0/3
% LIST 3: Ascending list of as-yet-untried doses
% de(<|=|>, [M|_], [V|Vs], [Z|Zs])
%% Hmmm... It seems I get some traction by treating LIST 2 as the concatenation
%% of LISTs 2 + 1 as described here, as if LIST 2 were a list *difference* Vs-Ms.
%% The constraint that Ms is a tail of Vs ought to be PROVEN.
%% What does de(C, Ms, [V|Vs], Xs) **MEAN**?
%% - length(Ms, MTD) would be a 'safe guess' for MTD, although possibly too low
%% - We have just done V~T with C giving a judgment via T as follows:
%%   (<) T was low enough to drive escalation
%%   (>) T was high enough to demand de-escalation
%%   (=) T was 'middling' value indicating another cohort at current dose V
%% - Xs are as-yet untried
% It seems to me that I am losing something of the CONTEXT of the C.
% Maybe this impression arises only because I have failed to clarify
% exactly WHEN the C applies! Is it a description of what has just
% occurred at the state (Ms, Vs, Xs)? YES!
% ---
% I think there may be some value to reordering the arguments to read
% in descending order, and reassigning the >,=,< conditionals to say
% which direction (greater, same, lower) we should SHIFT the current dose.
% But WHAT IS THE CURRENT DOSE in this predicate? Can I say it's the head
% of append(Vs, Ms)?
% What if I say the whole dose series is always append(Zs, Vs, Ms), and
% that the commas here are 'bookmarks' within the list?
current_dose(Vs, Ms, D) :-
    append(Vs, Ms, Ds),
    length(Ds, D).
de(>, [], Vs, Ms) --> [mtd_notfound(Top)],
		      { current_dose(Vs, Ms, Top) }. % 'hit ceiling'
de(>, [X|Xs], Vs, Ms) --> [D^T], % escalating to dose X ==> enrolling 3, observing T DLTs,
			 { current_dose([X|Vs], Ms, D),
			   tox(T),
			   zcompare(C, 1, T) },
			 de(C, Xs, [X|Vs], Ms).
de(=, Zs, Vs, Ms) --> [D-T], % NB: de(=,...) is a further cohort
		      { current_dose(Vs, Ms, D),
			tox(T) },
		      %zcompare(C, 1, T),
		      (	  {0 #= T}, de(>, Zs, [], Vs)
		      ;	  {0 #< T}, de(<, Vs, Ms)
		      ).
% Note the 'loss of arity' of de(<, ., .) ~~ "all downhill from here.."
de(<, [], Ms) --> [declare(MTD)], { length(Ms, MTD) }.
de(<, [V|Vs], Ms) --> [D:T],
		      { length([V|Vs],D),
			tox(T) },
		      (	  {T #= 0}, [declare(D)]
		      ;	  {T #> 0}, de(>, Ms, Vs)
		      ).

?- phrase(detrial(2), Tr).
%@ Tr = [1^0, 2^0, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^1, 2-0, mtd_notfound(2)] ;
%@ Tr = [1^0, 2^1, 2-_10184, 2:0, declare(2)],
%@ _10184 in 1..3 ;
%@ Tr = [1^1, 1-0, 2^0, mtd_notfound(2)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-0, mtd_notfound(1)] ;
%@ Tr = [1^1, 1-0, 2^1, 2-_14258, 1:0, declare(1)],
%@ _14258 in 1..3 ;
%@ Tr = [1^1, 1-_15762, 1:0, declare(1)],
%@ _15762 in 1..3 ;
%@ false.

/* Consider doing something to obtain a fair enumeration
 * by ordering the trials in descending order of toxicity.
 */

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

path_matrix_([declare_mtd(_)], _).
path_matrix_([mtd_notfound(_)], _).

path_matrix_([], (_,_)).

ground_or_nil(Term) :- ground(Term).
ground_or_nil('NA') :- true.

rep(E, N, L) :-
    (	N #= 0 -> L = []
    ;	N #> 0 -> (N_1 #= N - 1, L = [E|Es], rep(E, N_1, Es))
    ).

?- rep(a, 10, L).
%@ L = [a, a, a, a, a, a, a, a, a|...].
?- rep(X, N, [a,a,a,a,a]).
%@ false.
?- rep(a, 5, [a,a,a,a,a]).
%@ true.

repe(_, 0, []).
repe(E, N, [E|Es]) :-
    N #> 0,
    N_1 #= N - 1,
    repe(E, N_1, Es).

?- repe(X, N, [a,a,a,a,a]).
%@ X = a,
%@ N = 5 ;
%@ false.

?- repe(E, N, Es).
%@ N = 0,
%@ Es = [] ;
%@ N = 1,
%@ Es = [E] ;
%@ N = 2,
%@ Es = [E, E] ;
%@ N = 3,
%@ Es = [E, E, E] ;
%@ N = 4,
%@ Es = [E, E, E, E] .
% Dang!

%% Write out separate R input files for D in 2..7
write_escalation_array(D) :-
    atom_concat('T', D, Filename),
    atom_concat(Filename, '.tab', File),
    open(File, write, OS),
    findall(P, phrase(esc(0, 0..D), P), Paths),
    length(Paths, Len),
    write(Len),
    rep(D, Len, Ds),
    maplist(path_matrix, Paths, Ds, Ms),
    write_escalation_matrices(Ms, OS).

write_escalation_matrices([(C1,C2)|Ms], OS) :-
    write_matrix_row(C1, OS),
    write_matrix_row(C2, OS),
    write_escalation_matrices(Ms, OS).

write_escalation_matrices([], OS) :- close(OS).

write_matrix_row([E|Es], OS) :-
    write(OS, E),
    (	Es = [] -> nl(OS)
    ;	write(OS, '\t'),
	write_matrix_row(Es, OS)
    ).

/*
Some attempts to be even more explicit about the reasoning
being done in the trial//1 DCG above, while also achieving
greater generality such as will be required for describing
(e.g.) Simon &al accelerated titration phase.   
*/

% Initially, and strictly for convenience, I will continue
% to model the doses as integers. But this should yield to
% a data representation (like Peano) that supports stronger
% proofs of correctness.

% What is a valid 'top dose' in a prespecified dose range?
topdose(TopDose) :- TopDose in 1..10, indomain(TopDose).
%?- topdose(D).
%@ D = 1 ;
%@ D = 2 ;
%@ D = 3 ;
%@ D = 4 ;
%@ D = 5 ;
%@ D = 6 ;
%@ D = 7 ;
%@ D = 8 ;
%@ D = 9 ;
%@ D = 10.
%% NB: Oncology dose-finding trials in practice rarely
%%     pre-specify more than 10 distinct doses.

% What is a conceivable set of doses we might be considering
% (at any given moment in a dose-finding trial) as potential
% MTDs?
candidate_mtds([], TopDose) :- topdose(TopDose).
candidate_mtds([Singleton], TopDose) :-
    topdose(TopDose),
    Singleton in 1..TopDose,
    indomain(Singleton).
%?- candidate_mtds(C, TD).
%@ C = [],
%@ TD = 1 ;
%@ C = [],
%@ TD = 2 ;
%@ C = [9],
% ...as expected...
%@ TD = 10 ;
%@ C = [10],
%@ TD = 10.
candidate_mtds([Lowest,Next|Higher], TopDose) :-
    Lowest in 1..TopDose, indomain(Lowest),
    Lowest #= Next - 1, % Z is overkill for this!
    candidate_mtds([Next|Higher], TopDose).
%?- candidate_mtds(C, 3).
%@ C = [] ;
%@ C = [1] ;
%@ C = [2] ;
%@ C = [3] ;
%@ C = [1, 2] ;
%@ C = [1, 2, 3] ;
%@ C = [2, 3] ;
%@ false.

% That's an embarrassingly excessive use of Z.
% The essential notion is that of *succession*,
% for which lists suffice. What about even using
% a LIST DIFFERENCE?
% If we were currently considering [D|Ds] - Es,
% Es is a list, then we'd be saying that D is the
% lowest of ALL doses (D..inf) being considered,
% while Es are the doses we regard as excessive,
% say 5..inf. In such a case, we are (currently)
% actively considering the doses D..4.
% I sense distinctly the spirit of operating on
% the tape of a Turing machine!

toolow_maybe(Low, Maybe) :-
    TopDose in 1..10, indomain(TopDose),
    % How do I get the list [1,2,...,TopDose]?
    append(Low, Maybe, 1..TopDose).


