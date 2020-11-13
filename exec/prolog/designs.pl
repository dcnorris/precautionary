% Attempt an enumeration of ALL dose-finding designs over 2 doses

/* DISCUSS: How to switch on Prolog system?
*/

:- use_module(library(clpfd))  % SWI
   ;use_module(library(clpz)). % Scryer

%:- use_module(library(clpfd)).  % SWI
%:- use_module(library(clpz)). % Scryer
:- use_module(library(pio)).
:- use_module(library(pairs)).

%% I begin with designs that enroll cohorts of 3,
%% up to a maximum of 6 patients per dose.
%% TODO: Consider whether cohort size of 3 might
%%       emerge as a consequence of more basic
%%       considerations!
denominator(N1 - N2) :-
    N1 in 0 \/ 3 \/ 6,
    N2 in 0 \/ 3 \/ 6.
%    member(N1, [0, 3, 6]), % TODO: N1 in 0 \/ 3 \/ 6
%    member(N2, [0, 3, 6]).

%?- denominator(X - Y).
%@ X in 0\/3\/6,
%@ Y in 0\/3\/6.

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
% in a sense!

%% COMPARE:
%?- ( X#= 3 ; X #= 2).
%@ X = 3 ;
%@ X = 2.
%% VS:
%?- ( X#= 3 #\/ X #= 2).
%@ X in 2..3.

%?- tuples_in([[N1]], [[0],[3],[6]]).
%@    clpz:(N1 in 0\/3\/6)
%@ ;  false.
%@ N1 in 0\/3\/6. # SWI

%?- denominator(X - Y).
%@ X in 0\/3\/6,
%@ Y in 0\/3\/6.

:- debug.
%@ true.
numerator_denominator(T1 - T2, N1 - N2) :-
    denominator(N1 - N2),
    T1 #>= 0, T1 #=< N1,
    T2 #>= 0, T2 #=< N2.
%    T1 in 0..N1,
%    T2 in 0..N2.

%?- numerator_denominator(T, N).
%@ T = _51576-_51578,
%@ N = _51594-_51596,
%@ _51576 in 0..6,
%@ _51594#>=_51576,
%@ _51594 in 0\/3\/6,
%@ _51578 in 0..6,
%@ _51596#>=_51578,
%@ _51596 in 0\/3\/6.

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
unacceptable_tally(Q1 -  Q2) :-
    (	step_too_cautious(Q1, Q2)
    ;	step_too_bold(Q1, Q2)
    ).

%% Seeing we tallied up 0 DLTs at dose 1 without ever having
%% tried dose 2.
step_too_cautious(0/6, 0/0).

step_too_bold(T/3, 3/3) :- T #> 0.
step_too_bold(T/6, 3/3) :- T #> 1.

%?- unacceptable_tally(UT).
%@ UT = 0/6-0/0 ;
%@ UT = _85346/3-3/3,
%@ _85346 in 1..sup ;
%@ UT = _86692/6-3/3,
%@ _86692 in 2..sup.

/*
From the domain-expert perspective, it is far more
natural to describe tallies that are UNacceptable.

TODO: Exploit features such as #\ and tuple_in to
      achieve such an expression.
*/

acceptable_tally(T1 - T2, N1 - N2) :-
    T2 #< 5, % Want never to see 5+ toxicities at higher dose,
    T1 #< 4, % and never even 4+ toxicities at a lower dose.
    zcompare(C, N1, 6),
    acceptable_tally_(C, T1/N1 - T2/N2).

% What is acceptable in case N1 #= 6?
acceptable_tally_(=, T1/_ - T2/N2) :-
    (	T1 #> 0  %% i.e., N2 = 0 ==> T1 > 0, meaning that
    ;	N2 #> 0  %% we didn't dawdle at a too-low dose.
    ).

% What is acceptable in case N1 #< 6?
acceptable_tally_(<, T1/N1 - T2/N2) :-
    (	N2 #= 0  %% EITHER we haven't enrolled dose 2 yet...
    ;	T1 #= 0, N1 #>= 3 %% OR dose 1 had 0/3 tox tally.
    ).

%% NB: Clause for the '>' case is purposely absent:
%%acceptable_tally(>, _) :- false.

%?- acceptable_tally(T, N).
%@ T = _89184-_89186,
%@ N = 6-_89204,
%@ _89184 in 1..3,
%@ _89186 in inf..4 ;
%@ T = _91188-_91190,
%@ N = 6-_91208,
%@ _91188 in inf..3,
%@ _91190 in inf..4,
%@ _91208 in 1..sup ;
%@ T = _93552-_93554,
%@ N = _93570-0,
%@ _93552 in inf..3,
%@ _93554 in inf..4,
%@ _93570 in inf..5 ;
%@ T = 0-_95756,
%@ N = _95772-_95774,
%@ _95756 in inf..4,
%@ _95772 in 3..5.
%% TODO: Examine the above for correctness.

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

%@ true.
allowable_step(S1 - S2) :-
    denominator(N1a - N2a),
    denominator(N1b - N2b),
    N1b #= N1a + S1,
    N2b #= N2a + S2,
    S12 #= S1 + S2,
    S1 #>= 0,
    S2 #>= 0,
    S12 #< 4. % let's say, 4+ DLTs in 1 step is unacceptable

%?- denominator(N1a - N2a), denominator(N1b - N2b), N1b #= N1a + S1, N2b #= N2a + S2, S12 #= S1 + S2, S12 #< 4, S1 #>= 0, S2 #>= 0.
%@ N1a in 0\/3\/6,
%@ N1a+S1#=N1b,
%@ S1 in 0..3,
%@ S1+S2#=S12,
%@ S2 in 0..3,
%@ N2a+S2#=N2b,
%@ N2a in 0\/3\/6,
%@ N2b in 0\/3\/6,
%@ S12 in 0..3,
%@ N1b in 0\/3\/6.

%?- allowable_step(A - B).
%@ A in 0..3,
%@ A+B#=_20888,
%@ _20916+A#=_20912,
%@ B in 0..3,
%@ _20964+B#=_20960,
%@ _20964 in 0\/3\/6,
%@ _20960 in 0\/3\/6,
%@ _20888 in 0..3,
%@ _20916 in 0\/3\/6,
%@ _20912 in 0\/3\/6.

%% "Residual program"
%% The answer is a transformation of the program!

step(T1a/N1a - T2a/N2a, T1b/N1b - T2b/N2b) :-
    denominator(N1a - N2a),
    denominator(N1b - N2b),
    S1 #= N1b - N1a,
    S2 #= N2b - N2a,
    allowable_step(S1 - S2),
    acceptable_tally(T1b/N1b - T2b/N2b),
    acceptable_tally(T1a/N1a - T2a/N2a),
    true.

%?- step(A, B).

/* - - - -
%% TODO: Read memoization in Markus's book
:- dynamic memo/1.
memo(Goal) :- ( Goal -> true ; Goal, assertz(Goal)).

%% Dangeous territory!
:- dynamic memo_/1.
memo(Goal) :-
( memo_(Goal) -> true
; once(Goal),
assertz(memo_(Goal))
).
*/

%% Big questions: what are we talking about at all?

%% We can ABSORB the unification and backtracking!
%% Compare with the call stack of Lisp, which is implicit.

