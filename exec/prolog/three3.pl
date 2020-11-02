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
	  *indomain(T). % Can easly switch off labeling here!

% Mnemonic: * is ^ that 'splatted' on dose ceiling.
esc(Hi, Lo..Hi) --> [Hi * T], { Lo #< Hi,
				tox(T) },
		    (   { T #=< 1 }, [mtd_notfound(Hi)]
		    ;	{ T #>= 2 }, des(Hi, Lo)
		    ).
esc(D, Lo..Hi) --> [D1 ^ T], { D1 #= D + 1,
			       D1 in Lo..Hi,
			       tox(T) },
		   (   {T #= 0}, esc(D1, Lo..Hi)
		   ;   {T #= 1}, sta(D1, Lo..Hi)
		   ;   {T #> 1}, des(D1, Lo)
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
	       (   { D_1 #= Lo }, [declare_mtd(Lo)]
	       ;   { D_1 #> Lo }, [D_1 : T], { tox(T) },
		   (   { T #=< 1 }, [declare_mtd(D_1)]
		   ;   { T #>= 2 }, des(D_1, Lo)
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
detrial(N) --> { length(Ds, N) }, de(<, [], [], Ds).
% de(C, D, E) means we just observed a tox count (<,=,>)~(low,middling,high)
% on dose length(D), with additional doses E above it still being considered.
% Perhaps the full computation requires maintaining 3 lists, the first of which
% is doses known to be safe, the second being those tried only once (and therefore
% not presumed safe or unsafe) and the third being doses not yet tried.
de(<, Ls, []) --> [mtd_notfound(Top)], { length(Ls, Top) }. % 'hit ceiling'
de(<, Ls, [R|Rs]) --> [D^T], % escalating to dose R ==> enrolling 3, observing T DLTs,
		      { length([R|Ls], D), T in 0..3 },
		      zcompare(C, T, 1),
		      de(C, [R|Ls], Rs).
de(=, [X|Ls], Rs) --> [D-T], % NB: de(=,...) is a further cohort
		      { length([X|Ls], D), T in 0..3 },
		      %%zcompare(C, T, 1),
		      (   T #< 1, de(<, [X|Ls], Rs)
		      ;   declare(T, X)
		      ).
de(>, [X|Ls], Rs) --> [D:T], { T in 2..3 }. % NB: de(>,...) is a 2nd cohort on 'downslope'
			


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


% I also need to consider the manner of representing DOSES themselves!
% This might well be a LIST of dose ALIQUOTS, such that the 'final dose'
% is constituted by the SUM over the list. Note how this representation
% (or 'interpretation'?) automatically builds in the monotonicity or
% sorting of the doses.

% At any point, the STATE of the trial is a pair of lists As ^ [A|More],
% where the current dose 'under consideration' is to be understood as
% sum(As) + A, and there are possibly additional higher doses that can
% be 'constructed' ('administered'!) by drawing up the further aliquots
% in the list More.
% To the extent it proves desirable to instantiate the solutions, I can
% very well begin with simple lists of the term 'a'!
%
% OTOH, maybe I don't even need the dose designations (as in 'a') apart
% from the position in the list! If I discard the idea of truncating the
% list as a means to indicate that some doses have been ruled out-of-bounds,
% then I can simply preserve all information about the doses tried and the
% DLT fractions observed. There is no need to 'chop off' the Turing tape,
% after all!

% What does a  cohort look like?
cohort(DLTs/N) :-
    N #= 3, % initially, we permit only cohorts of 3
    DLTs in 0..N, indomain(DLTs).
%?- cohort(C).
%@ C = 0/3 ;
%@ C = 1/3 ;
%@ C = 2/3 ;
%@ C = 3/3.

% What does a DLT tally look like?
tally(DLTs/Enrolled) :-
    Enrolled in 0..6, indomain(Enrolled),
    % -------------------------------------------------------------------------------
    Ncohorts in 0..2,         % TODO: Release this constraint to permit finer-grained
    indomain(Ncohorts),       %       enrollment 0..6 as above, which will be needed
    Enrolled #= Ncohorts * 3, %       to model Simon &al 1997 accelerated titration.
    % -------------------------------------------------------------------------------
    DLTs in 0..Enrolled, indomain(DLTs).
%?- tally(C).
%@ C = 0/0 ;
%@ C = 0/3 ;
%@ C = 1/3 ;
%@ C = 2/3 ;
%@ C = 3/3 ;
%@ C = 0/6 ;
%@ C = 1/6 ;
%@ C = 2/6 ;
%@ C = 3/6 ;
%@ C = 4/6 ; % Note that none of these 3 tallies
%@ C = 5/6 ; % could ever arise in a 3+3 trial.
%@ C = 6/6.  % Not asking this predicate to do it all!

enrollable_tally(T0/N0) :-
    tally(T0/N0),
    (	N0 #= 0
    ;	N0 #= 3, T0 #=< 1
    ).
%?- enrollable(Q).
%@ Q = 0/0 ;
%@ Q = 0/3 ;
%@ Q = 1/3 ;
%@ false.

% How do cohorts ACCUMULATE?
tally0_cohort_tally(T0/N0, T_/N_, T/N) :-
    enrollable_tally(T0/N0),
    cohort(T_/N_),
    tally(T/N),
    T #= T0 + T_,
    N #= N0 + N_.
%?- tally0_cohort_tally(T0, C, T).
%@ T0 = 0/0,
%@ C = T, T = 0/3 ;
%@ T0 = 0/0,
%@ C = T, T = 1/3 ;
%@ T0 = 0/0,
%@ C = T, T = 2/3 ;
%@ T0 = 0/0,
%@ C = T, T = 3/3 ;
%@ T0 = C, C = 0/3,
%@ T = 0/6 ;
%@ T0 = 0/3,
%@ C = 1/3,
%@ T = 1/6 ;
%@ T0 = 0/3,
%@ C = 2/3,
%@ T = 2/6 ;
%@ T0 = 0/3,
%@ C = 3/3,
%@ T = 3/6 ;
%@ T0 = 1/3,
%@ C = 0/3,
%@ T = 1/6 ;
%@ T0 = C, C = 1/3,
%@ T = 2/6 ;
%@ T0 = 1/3,
%@ C = 2/3,
%@ T = 3/6 ;
%@ T0 = 1/3,
%@ C = 3/3,
%@ T = 4/6 ;
%@ false.

% How do tallies COMPARE?
safer_than(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0, % NB: These N>0 conditions are implicit in the #< below,
    tally(T1/N1), N1 #> 0, %     but are stated explicitly here to expose the logic.
    T0*N1 #< N0*T1. % 'cross-multiply'
%?- safer_than(Q0, Q1).
%@ Q0 = 0/3,
%@ Q1 = 1/3 ;
% ...  as expected ...
%@ Q0 = 5/6,
%@ Q1 = 3/3 ;
%@ Q0 = 5/6,
%@ Q1 = 6/6 ;
%@ false.
%?- safer_than(Q, 1/3).
%@ Q = 0/3 ;
%@ Q = 0/6 ;
%@ Q = 1/6 ;
%@ false.

noworse_than(T0/N0, T1/N1) :-
    tally(T0/N0), N0 #> 0, % NB: By contrast with the N>0 conditions in safer_than/2,
    tally(T1/N1), N1 #> 0, %     these are strictly necessary for correctness here
    T0*N1 #=< N0*T1.       % <-- because the weaker comparison '#=<' is used.
%?- noworse_than(Q0, Q1), \+ safer_than(Q0, Q1), \+ Q0 = Q1. % find nontrivial set difference
%@ Q0 = 0/3,
%@ Q1 = 0/6 ;
%@ Q0 = 1/3,
%@ Q1 = 2/6 ;
%@ Q0 = 2/3,
%@ Q1 = 4/6 ;
%@ Q0 = 3/3,
%@ Q1 = 6/6 ;
%@ Q0 = 0/6,
%@ Q1 = 0/3 ;
%@ Q0 = 2/6,
%@ Q1 = 1/3 ;
%@ Q0 = 4/6,
%@ Q1 = 2/3 ;
%@ Q0 = 6/6,
%@ Q1 = 3/3 ;
%@ false.
%?- noworse_than(Q, 1/6).
%@ Q = 0/3 ;
%@ Q = 0/6 ;
%@ Q = 1/6 ;
%@ false.

:- op(900, xfx, user:(&=<)).
&=<(Q1, Q2) :- noworse_than(Q1, Q2).
%?- 1/3 &=< Q.
%@ Q = 1/3 ;
%@ Q = 2/3 ;
%@ Q = 3/3 ;
%@ Q = 2/6 ;
%@ Q = 3/6 ;
%@ Q = 4/6 ;
%@ Q = 5/6 ;
%@ Q = 6/6 ;
%@ false.

:- op(900, xfx, user:(&<)).
&<(Q1, Q2) :- safer_than(Q1, Q2).
%?- Q &< 1/3.
%@ Q = 0/3 ;
%@ Q = 0/6 ;
%@ Q = 1/6 ;
%@ false.

% Describe a LIST of TALLIES for consecutive prespecified doses:
tallies([]).
tallies([Q|Qs]) :-
    length(Qs, _), % enumerate solutions fairly
    tally(Q),
    tallies(Qs).
%?- length(C, 1), tallies(C).
%@ C = [0/0] ;
%@ C = [0/3] ;
%@ C = [1/3] ;
%@ C = [2/3] ;
%@ C = [3/3] ;
%@ C = [0/6] ;
%@ C = [1/6] ;
%@ C = [2/6] ;
%@ C = [3/6] ;
%@ C = [4/6] ; % Note that none of these cases
%@ C = [5/6] ; % could ever happen in a 3+3 trial.
%@ C = [6/6].  % This predicate can't do it all!
%?- length(C, 2), tallies(C).
%@ C = [0/0, 0/0] ;
%@ C = [0/0, 0/3] ;
%@ C = [0/0, 1/3] ;
%@ C = [0/0, 2/3] ;
% ... as expected ...
%@ C = [6/6, 4/6] ;
%@ C = [6/6, 5/6] ;
%@ C = [6/6, 6/6]. 
%?- length(C, 2), tallies(C), false. % Does it terminate?
%@ false.                            % Yes.

% Convenient access to the list heads would seem to advise holding the
% current state of the trial as a PAIR of lists, the first descending
% and the second ascending. That is, they are two STACKS side-by-side,
% with the top (front) elements being (on the left) the highest of the
% lower part of doses, and (on the right) the lowest among the higher
% doses.
% Let me also aim to separate the stack-shuffling ACTIONS from the
% subsequent (automatic) ENROLLMENT vs STOPPING of the trial.
% After each ACTION, the trial is ready for CONTINUATION = ENROLL | STOP.
% If ENROLLMENT is to take place, then it is at the top of the right stack.
% In order to separate those events, I will encode the alternating state
% of the trial using ^ and : infix operators. The mnemonic content is that
% the ^ looks like the scales of Justice, with some kind of 'weighing' or
% judgment made about the doses, while : looks like '...' indicating the
% pending continuation of the trial.
% ***
% I'm doubling down now on my view that (escalate | stay | deescalate)
% register MAINTAINANCE OF OUR ATTENTION on the most interesting dose.
% This is separate from what we DECIDE TO DO after inspecting this
% currently-most-interesting dose. In theory, I could even rename half
% of these functors to state0_FOCUS_state().
state0_action_state(Ls ^ [Q|Rs], escalate, [Q|Ls] : Rs) :-
    Q &< 1/3.
state0_action_state(Ls^[Q|Rs], stay, Ls:[Q|Rs]) :-
    Q = 1/3.
% TODO: Express this clause more ABSTRACTLY, then later PROVE
%       that deescalation occurs only to a pre-existing 0/3 tally.
state0_action_state([0/3|Ls] ^ Rs, deescalate, Ls : [0/3|Rs]) :-
    Rs = [Q|_],
    1/3 &< Q.
/****
The following general queries seem to PROVE
key properties of 3+3 have been attained.
 ****/
%?- state0_action_state(S0, deescalate, S).
%@ S0 = [0/3|_47804]^[2/3|_47840],
%@ S = [_47804]:[0/3, 2/3|_47840] ;
%@ S0 = [0/3|_47804]^[3/3|_47840],
%@ S = [_47804]:[0/3, 3/3|_47840] ;
%@ S0 = [0/3|_47804]^[3/6|_47840],
%@ S = [_47804]:[0/3, 3/6|_47840] ;
%@ S0 = [0/3|_47804]^[4/6|_47840],
%@ S = [_47804]:[0/3, 4/6|_47840] ;
%@ S0 = [0/3|_47804]^[5/6|_47840],
%@ S = [_47804]:[0/3, 5/6|_47840] ;
%@ S0 = [0/3|_47804]^[6/6|_47840],
%@ S = [_47804]:[0/3, 6/6|_47840] ;
%@ false.
%?- state0_action_state(S0, stay, S).
%@ S0 = _61938^[1/3|_61946],
%@ S = [_61938]:[1/3|_61946].
%?- state0_action_state(S0, escalate, S).
%@ S0 = _64680^[0/3|_64688],
%@ S = [0/3|_64680]:_64688 ;
%@ S0 = _64680^[0/6|_64688],
%@ S = [0/6|_64680]:_64688 ;
%@ S0 = _64680^[1/6|_64688],
%@ S = [1/6|_64680]:_64688 ;
%@ false.

% NOTE: This currently quite trivial enrollment step may well grow
%       more complex with introduction of an accelerated titration
%       phase or other such variations on standard 3+3.
state0_action_state(Ls : [Q0 | Rs], enroll, Ls ^ [Q | Rs]) :-
    tally0_cohort_tally(Q0, _, Q).
%?- state0_action_state(S0, enroll, S).
%@ S0 = _4386:[0/0|_4394],
%@ S = _4386^[0/3|_4394] ;
%@ S0 = _4386:[0/0|_4394],
%@ S = _4386^[1/3|_4394] ;
%@ S0 = _4386:[0/0|_4394],
%@ S = _4386^[2/3|_4394] ;
%@ S0 = _4386:[0/0|_4394],
%@ S = _4386^[3/3|_4394] ;
%@ S0 = _4386:[0/3|_4394], % NB: In a 3+3 trial, _:[0/3|_] shouldn't
%@ S = _4386^[0/6|_4394] ; %     ever appear as a state0. But I think
%@ S0 = _4386:[0/3|_4394], %     it's NOT this predicate's proper role
%@ S = _4386^[1/6|_4394] ; %     to impose this condition.
%@ S0 = _4386:[0/3|_4394],
%@ S = _4386^[2/6|_4394] ;
%@ S0 = _4386:[0/3|_4394],
%@ S = _4386^[3/6|_4394] ;
%@ S0 = _4386:[1/3|_4394],
%@ S = _4386^[1/6|_4394] ;
%@ S0 = _4386:[1/3|_4394],
%@ S = _4386^[2/6|_4394] ;
%@ S0 = _50:[1/3|_58],
%@ S = _50^[3/6|_58] ;
%@ S0 = _50:[1/3|_58],
%@ S = _50^[4/6|_58] ;
%@ false.

% If we 'take seriously' the analogy with stacks, then you don't
% realize the right-hand stack had no further doses until you pop
% the top element off and notice 'nothing there'. On that account,
% it seems natural to make the determination "I *CAN'T* continue!"
% on the ':' side.
% (The sheer simplicity of this predicate also argues in favor.)
% (Can I appeal also to the 'exceptional' nature of this determination?)
state0_action_state(Ls : [], stop, mtd_notfound(D)) :-
    length(Ls, D). % RP2D is highest prespecified dose, D.

/*
Instead of regarding merely TALLIES as presumably safe/toxic,
we apply these notions rather (and more properly) to the DOSES
themselves -- as represented through lists.
Interestingly, this helps render state0_action_state( _ , stop, declare_mtd)
through one goal without a special case.
*/
presumably_safe([]). % [] corresponds to zero dose
presumably_safe([Q|_]) :- % "no need to enroll further at this dose to say it's okay"
    tally(Q),
    Q = T/6,   % TODO: Might I obtain the no-deescalation variant of 3+3
    T in 0..1. %       simply by modifying these 2 goals?

presumably_toxic([Q|_]) :-
    tally(Q),
    Q = T/N,
    T #> 1.

%% About which (indeterminate) tallies do we presume nothing?
%?- tally(I), \+ presumably_safe([I]), \+ presumably_toxic([I]).
%@ I = 0/0 ;
%@ I = 0/3 ;
%@ I = 1/3 ;
%@ false.

state0_action_state(Ls : Rs, stop, declare_mtd(MTD)) :-
    presumably_safe(Ls),
    presumably_toxic(Rs),
    length(Ls, MTD).

actions(_) --> [].
actions(S0) --> [A->S],
        { state0_action_state(S0, A, S) },
        actions(S).

%% Let's examine conduct of a trial with 2 dose levels.
%% (Trials always start off like [0/0, ..., 0/0] : [].)
%?- phrase(actions([] : [0/0, 0/0]), Trial).
%@ Trial = [] ;
%@ Trial = [(enroll->[]^[0/3, 0/0])] ;
%@ Trial = [(enroll->[]^[0/3, 0/0]),  (escalate->[0/3]:[0/0])] ;
%@ Trial = [(enroll->[]^[0/3, 0/0]),  (escalate->[0/3]:[0/0]),  (enroll->[0/3]^[0/3])] ;
%@ Trial = [(enroll->[]^[0/3, 0/0]),  (escalate->[0/3]:[0/0]),  (enroll->[0/3]^[0/3]),  (escalate->[0/3, 0/3]:[])] ;
%@ Trial = [(enroll->[]^[0/3, 0/0]),  (escalate->[0/3]:[0/0]),  (enroll->[0/3]^[0/3]),  (escalate->[0/3, 0/3]:[]),  (stop->mtd_notfound(2))] .
%
%% FASCINATING! I had never thought about this case, but the decision of the program
%% in retrospect looks reasonable. I do rather suspect that the spirit of the '=< 1/6'
%% definition would require that any dose 'reported out' of the trial (e.g. as RP2D)
%% ought to have been tested in 6 patients.

/* TODO ...

1. Let the 3+3 'rules' be DEMONSTRATED by complete *singleton* solutions
   to a sequence of Prolog QUERIES. This is the 'proof of correctness'
   that I was supposing I had 'approximately' obtained with the matching
   solution-set sizes for trial//1 and esc//2.

2. The RATIONALITY of any 3+3 program will be tested by how readily it
   supports generalization to other variants and extensions.

a) Demonstrate that the alternate 3+3 design (without de-escalation) can
   readily be obtained by a small (ideally, *parametrized*) modification.

b) Similarly introduce the accelerated titration. Ideally here, we might
   demonstrate COMPOSITIONALITY such that the 3+3 with an A.T. phase can
   be expressed by sequencing within a DCG.

*/
