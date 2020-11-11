% Attempt an enumeration of ALL dose-finding designs over 2 doses
:- use_module(library(clpz)).
%@ caught: error(existence_error(source_sink,library(clpfd)),use_module/1)
:- use_module(library(pio)).
%@    true.
:- use_module(library(pairs)).
%@    true.
%@ true.

/* - - - - - 

*/

%% I begin with designs that enroll cohorts of 3,
%% up to a maximum of 6 patients per dose.
%% TODO: Consider whether cohort size of 3 might
%%       emerge as a consequence of more basic
%%       considerations!
denominator((N1, N2)) :-
    member(N1, [0, 3, 6]), % TODO: N1 in 0 \/ 3 \/ 6
    member(N2, [0, 3, 6]).
%?- #\ X in 1..5.
%@ X in inf..0\/6..sup.

%% TODO: Try negating tuples_in using #\.
%%       But be careful, as this is a less-used
%%       part of clpfd.

%% Interaction between tabling and constraints here
%% will probably be safe here, but might not be so
%% effective.

% You don't get constraints you can reason about TOGETHER
% when you use member/2 -- which 'hard-codes' the search,
% ib a sense!
% This is the connection between constraints and regular
% prolog programming.

/*
COMPARE:
( X#= 3 ; X #= 2)
( X#= 3 #\/ X #= 2)
*/

%?- X in 0..3, member(X, [0,2]).
%@ X = 0 ;
%@ X = 2.

%?- tuples_in([N1], [[0],[3],[6]]).
%@ ERROR: Arguments are not sufficiently instantiated
%@ ERROR: In:
%@ ERROR:   [18] throw(error(instantiation_error,_30188))
%@ ERROR:   [10] clpfd:tuples_in([_30224],[[0],...|...]) at /usr/local/Cellar/swi-prolog/8.2.1/libexec/lib/swipl/library/clp/clpfd.pl:4063
%@ ERROR:    [9] <user>
%@ ERROR: 
%@ ERROR: Note: some frames are missing due to last-call optimization.
%@ ERROR: Re-run your program in debug mode (:- debug.) to get more detail.
%@ ERROR: Arguments are not sufficiently instantiated
%@ ERROR: In:
%@ ERROR:   [18] throw(error(instantiation_error,_26636))
%@ ERROR:   [10] clpfd:tuples_in([_26672],[[0],...|...]) at /usr/local/Cellar/swi-prolog/8.2.1/libexec/lib/swipl/library/clp/clpfd.pl:4063
%@ ERROR:    [9] <user>
%@ ERROR: 
%@ ERROR: Note: some frames are missing due to last-call optimization.
%@ ERROR: Re-run your program in debug mode (:- debug.) to get more detail.

%% NONDETERMINISTIC because of member/2
%?- denominator((X,Y)).
%@ X = Y, Y = 0 ;
%@ X = 0,
%@ Y = 3 ;
%@ X = 0,
%@ Y = 6 ;
%@ X = 3,
%@ Y = 0 ;
%@ X = Y, Y = 3 ;
%@ X = 3,
%@ Y = 6 ;
%@ X = 6,
%@ Y = 0 ;
%@ X = 6,
%@ Y = 3 ;
%@ X = Y, Y = 6.

numerator_denominator((T1,T2), (N1,N2)) :-
    denominator((N1,N2)),
    T1 in 0..N1,
    T2 in 0..N2.

%?- numerator_denominator(T, N).
%@ T = N, N =  (0, 0) ;
%@ T =  (0, _30674),
%@ N =  (0, 3),
%@ _30674 in 0..3 ;
%@ T =  (0, _32320),
%@ N =  (0, 6),
%@ _32320 in 0..6 ;
%@ T =  (_33982, 0),
%@ N =  (3, 0),
%@ _33982 in 0..3 ;
%@ T =  (_35726, _35728),
%@ N =  (3, 3),
%@ _35726 in 0..3,
%@ _35728 in 0..3 ;
%@ T =  (_37658, _37660),
%@ N =  (3, 6),
%@ _37658 in 0..3,
%@ _37660 in 0..6 ;
%@ T =  (_39510, 0),
%@ N =  (6, 0),
%@ _39510 in 0..6 ;
%@ T =  (_41254, _41256),
%@ N =  (6, 3),
%@ _41254 in 0..6,
%@ _41256 in 0..3 ;
%@ T =  (_43186, _43188),
%@ N =  (6, 6),
%@ _43186 in 0..6,
%@ _43188 in 0..6. 
% Good.

/* - - -
Note that numerator_denominator/2 enumerates all possible toxicity tallies
in a 2-dose design, without regard for any characteristics of that design
that would render it suitable for any application whatsoever!
*/


/* - - - 

We end up with a DENOMINATOR LATTICE that looks like this:

6  _.._.._
   .......
   .......
3  _.._.._
   .......
   .......
0  _.._.._
   0  3  6

The explicitly-imposed cohort-of-3 constraint means that only
the 3*3=9 underscored points in the lattice are 'in play'.

At any point in this denominator lattice, there will be a feasible
set of tallies. Some of these tallies will however be UNACCEPTABLE,
meaning that we disallow any design under which they are reachable.
WLOG, then, we may confine our attention to the sets of ACCEPTABLE
TALLIES 'located' at each point of the lattice according to their
denominators.

Incidentally, it is worth noting that rotating the above lattice
clockwise by 135 degrees yields a TREE in which leftward movement
corresponds to choosing the lower dose and rightward branching to
choosing the higher dose. (Also, observe that these are not the
same as 'escalation' or 'de-escalation' because the latter concepts
are PATH-DEPENDENT.)

*/

/*
Certain tallies may be 'unacceptable', in a medicolegal sense.
More generally the principle is one of extreme 'regret', but
I do think the 'defensive' posture of medicolegalism expresses
the emotions at work here.

NB: Unlike the 'tally' of aliquots.pl, which corresponded to
a quotient T/N, here a 'tally' is the quotient PAIR for both
doses. 
*/
unacceptable_tally((Q1, Q2)) :-
    (	step_too_cautious(Q1, Q2)
    ;	step_too_bold(Q1, Q2)
    ).

%% Seeing we tallied up 0 DLTs at dose 1 without ever having
%% tried dose 2.
step_too_cautious(0/6, 0/0).

step_too_bold(T/3, 3/3) :- T #> 0.
step_too_bold(T/6, 3/3) :- T #> 1.

%?- unacceptable_tally(UT).
%@ UT =  (0/6, 0/0) ;
%@ UT =  (_20518/3, 3/3),
%@ _20518 in 1..sup ;
%@ UT =  (_21714/6, 3/3),
%@ _21714 in 2..sup.

/*
DISCUSS: From the domain-expert perspective, it is far more
natural to describe tallies that are UNacceptable. Are there
techniques or idioms (beyond zcompare/3) that facilitate a
'subtractive' style of expression of CLPZ-type constraints
without recourse to impure constructs ->/2 or \+? Would such
requirements imply developing a domain-specific language?
tor((T1,T2), (N1,N2)),
    T2 #< 5, % Want never to see 5+ toxicities at higher dose,
    T1 #< 4, % and never even 4+ toxicities at a lower dose.
    zcompare(C, N1, 6),
    acceptable_tally_(C, (T1/N1, T2/N2)).

% What is acceptable in case N1 #= 6?
acceptable_tally_(=, (T1/_, T2/N2)) :-
    (	T1 #> 0  %% i.e., N2 = 0 ==> T1 > 0, meaning that
    ;	N2 #> 0  %% we didn't dawdle at a too-low dose.
    ).

% What is acceptable in case N1 #< 6?
acceptable_tally_(<, (T1/N1, T2/N2)) :-
    (	N2 #= 0  %% EITHER we haven't enrolled dose 2 yet...
    ;	T1 #= 0, N1 #>= 3 %% OR dose 1 had 0/3 tox tally.
    ).

%% NB: Clause for the '>' case is purposely absent:
%%acceptable_tally(>, _) :- false.

/*

A DESIGN then consists of a set of arrows connecting TALLIES to
TALLY SETS that are reachable with a single, allowable step.
Provided that 'acceptability' constitutes a PREFERENCE RELATION
that is monotonic in toxicity (more toxicity is always worse),
then I think we may WLOG draw the arrows simply to the WORST
POSSIBLE OUTCOME TALLY.

Visualized as overlaid on the denominator lattice, then, the ARROWS
OF A DESIGN map TALLIES --> TALLIES, indicating for each escalation
decision the worst possible outcome. This is a STRONG FORMULATION
that FROM THE OUTSET excludes arrows corresponding even to very rare
catastrophic events. But this is just fine, since formulating the
problem this way expresses OUR INTENT TO CONDUCT TRUE SAFETY ANALYSIS.

Of note, while we do obtain a CATEGORY here, in which the tallies
(T1/N1, T2/N2) are the OBJECTS and the steps are the ARROWS, unless
category theory has some notion of 'primitive' arrows then we lose
the important distinction between (a) multi-step arrows formed by
the composition of the steps---during which information emerges and
decisions are made---and (b) single 'atomic' steps. This INFORMATIONAL
aspect is all-important!

*/

:- table allowable_step/1.
%@ true.
allowable_step((S1,S2)) :-
    S1 #>= 0,
    S2 #>= 0,
    denominator((N1a,N2a)),
    denominator((N1b,N2b)),
    N1b #= N1a + S1,
    N2b #= N2a + S2,
    S1 + S2 #< 4. % let's say, 4 DLTs in 1 step is unacceptable

%?- allowable_step(S).
%@ S =  (3, 0) ;
%@ S =  (0, 3) ;
%@ S =  (0, 0).

step((T1a/N1a,T2a/N2a), (T1b/N1b,T2b/N2b)) :-
    denominator((N1a, N2a)),
    denominator((N1b, N2b)),
    S1 #= N1b - N1a,
    S2 #= N2b - N2a,
    allowable_step((S1,S2)),
    acceptable_tally((T1b/N1b,T2b/N2b)),
    acceptable_tally((T1a/N1a,T2a/N2a)),
    true.

