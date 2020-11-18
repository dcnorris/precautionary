% Attempt an enumeration of ALL dose-finding designs over 2 doses

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

/* I have moved these to implementation-specific init files
:- use_module(library(clpfd)).  % SWI
:- use_module(library(clpz)).   % Scryer
*/

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
%@ N1 in 0\/3\/6.

%?- denominator(X - Y).
%@ X in 0\/3\/6,
%@ Y in 0\/3\/6.

:- debug.
%@ true.
numerator_denominator(T1 - T2, N1 - N2) :-
    denominator(N1 - N2),
    T1 #>= 0, T1 #=< N1,
    T2 #>= 0, T2 #=< N2.

%?- numerator_denominator(T, N).
%@ T = _39444-_39446,
%@ N = _39462-_39464,
%@ _39444 in 0..6,
%@ _39462#>=_39444,
%@ _39462 in 0\/3\/6,
%@ _39446 in 0..6,
%@ _39464#>=_39446,
%@ _39464 in 0\/3\/6.

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
corresponds to choosing the lower dose for the next cohort, and
rightward branching to choosing the higher dose for the next cohort.
(Also, observe that these are DISTINCT FROM 'escalation' and
'de-escalation' because the latter concepts are PATH-DEPENDENT.)

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
%@ UT = _42196/3-3/3,
%@ _42196 in 1..sup ;
%@ UT = _43542/6-3/3,
%@ _43542 in 2..sup.

/*
From the domain-expert perspective, it is far more
natural to describe tallies that are UNacceptable.

Accordingly, we use #\ in the following.

TODO: Remain alert for opportunities to use tuples_in/2,
      especially with the generalization to D > 2.

TODO: Consider whether describing a partial order
      on the tallies (worse-than) would help with
      the formulation of 'acceptable'.
*/

%% TODO: Consider 

acceptable_tally(T1/N1 - T2/N2) :-
    numerator_denominator(T1 - T2, N1 - N2),
    % TODO: Articulate a rationale for these 2 'never' demands:
    #\ T2 in 5..sup, % Want never to see 5+ toxicities
    #\ T1 in 5..sup, % at either dose.
    % Now, depending on net enrollment (so far) at the lower dose,
    % we have some further things to say about what is acceptable...
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

%?- acceptable_tally(T1/N1 - T2/N2).
%@ N1 = 6,
%@ T1 in 1..4,
%@ T2 in 0..4,
%@ N2#>=T2,
%@ N2 in 0\/3\/6 ;
% Acceptable worst-case tallies under this solution are:
% 1. 4/6 - 4/6
% 2. 4/6 - 3/3
% 3. 4/6 - 0/0
% These are ALL CORRECT.
% At the other extreme, we have:
% 4. 1/6 - 0/* , which are fine too.
%@ N1 = 6,
%@ T1 in 0..4,
%@ T2 in 0..4,
%@ N2#>=T2,
%@ N2 in 3\/6 ;
% Relative to the solution above, this solution appropriately
% disallows 0/6 - 0/0, which is indeed undesirable as it sees
% a full 0/6 toxicities at low dose without having escalated.
%@ T2 = N2, N2 = 0,
%@ T1 in 0..3,
%@ N1#>=T1,
%@ N1 in 0\/3 ;
% What these cases allow is 0/0 - 0/0 (okay) and */3 - 0/0.
% These are both fine, as they represent 'just starting'.
%@ T1 = 0,
%@ N1 = 3,
%@ T2 in 0..4,
%@ N2#>=T2,
%@ N2 in 0\/3\/6.
% These (finally) are the cases where 0/3 leads to escalation:
% 0/3 - 0/0 or */3 or {0,1}/3+*/3.
%% DONE: Examine the above for correctness.

/*

A DESIGN then consists of a set of arrows connecting TALLIES to
TALLY SETS that are reachable with a single, allowable step.
               vvvvv WRONG IDEA vvvvv
Provided that 'acceptability' constitutes a PREFERENCE RELATION
that is monotonic in toxicity (more toxicity is always worse),
then I think we may WLOG draw the arrows simply to the WORST
POSSIBLE OUTCOME TALLY.
               ^^^^^ WRONG IDEA ^^^^^

WAIT! Acceptability MAY VERY WELL NOT be fully monotonic in toxicity.
For example, some very low tallies may indicate a failure to escalate
fast enough. Accordingly, I require perhaps a more general notion of
what a design is. Yes, this is true!

So in fact I see CLP must ACTUALLY be brought into play SUBSTANTIVELY.
This is a GOOD THING! But it does mean that my formulation remains 
somewhat complex. The targets of arrows are perhaps mere denominators,
but the CONDITIONS on VALID arrows are UNIVERSALLY QUANTIFIED.

That is, an arrow is valid iff ALL POSSIBLE OUTCOME TALLIES are acceptable.

Visualized as overlaid on the denominator lattice, then, the ARROWS
OF A DESIGN map TALLIES --> DENOMINATORS, such that the denominator
*difference* indicates for each escalation decision how many to enroll
at each dose. (I am leaving open the possibility of mixed enrollment,
to retain generality.) The valid arrows are those for which every possible
end-tally is acceptable. Tallies that are reachable from a valid arrow,
but which form the tail of no valid arrow, are terminal states of the
design.

If I wanted to obtain a true CATEGORY, I might do this by recognizing
that in general the OBJECTS would be SETS OF TALLIES.

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

%?- allowable_step(A - B).
%@ A in 0..3,
%@ A+B#=_12040,
%@ _12068+A#=_12064,
%@ B in 0..3,
%@ _12116+B#=_12112,
%@ _12116 in 0\/3\/6,
%@ _12112 in 0\/3\/6,
%@ _12040 in 0..3,
%@ _12068 in 0\/3\/6,
%@ _12064 in 0\/3\/6.
% Interesting to appreciate that the unlabelled
% residual program isn't too helpful for vetting
% in this particular instance.
%?- allowable_step(A - B), indomain(A), indomain(B).
%@ A = B, B = 0 ; %% This is okay, if the identity arrow is useful.
%@ A = 0,  %% This is the case of
%@ B = 3 ; %% enrolling dose-level 2.
%@ A = 3,  %% Enrolling
%@ B = 0.  %% dose-level 1.

%?- X in 0..1, Y in 0..1, #\(X #= Y, Y #= 1), indomain(X), indomain(Y).
%@ X = Y, Y = 0 ;
%@ X = 0,
%@ Y = 1 ;
%@ false.

%?- X in 0..1, Y in 0..1, X + Y #< 2, indomain(X), indomain(Y).
%@ X = Y, Y = 0 ;
%@ X = 0,
%@ Y = 1 ;
%@ X = 1,
%@ Y = 0.

%% Here, we detail CONCEIVABLE OUTCOMES based on ARITHMETIC.
tally0_tally(T1a/N1a - T2a/N2a, T1b/N1b - T2b/N2b) :-
    acceptable_tally(T1a/N1a - T2a/N2a), % No sense discussing an unacceptable tally0, tho.
    denominator(N1b - N2b),
    E1 #= N1b - N1a,
    E2 #= N2b - N2a,
    E1 #>= 0, % Enrollment must be
    E2 #>= 0, % non-negative!
    E1 + E2 #> 0, % Net enrollment must be nontrivial.
    indomain(E1), % Avoid 'arguments not sufficiently
    indomain(E2), % instantiated' error
    T1new in 0..E1, % TODO: Try plainer constraints here,
    T2new in 0..E2, % to avoid premature labeling above.
    T1b #= T1a + T1new,
    T2b #= T2a + T2new.

%?- tally0_tally(0/3 - 0/0, T1/N1 - T2/N2), indomain(N1), indomain(N2).
%@ T1 = 0,
%@ N1 = N2, N2 = 3,
%@ T2 in 0..3 ;
%@ T1 = 0,
%@ N1 = 3,
%@ N2 = 6,
%@ T2 in 0..6 ;
%@ N1 = 6,
%@ T2 = N2, N2 = 0,
%@ T1 in 0..3 ;
%@ N1 = 6,
%@ N2 = 3,
%@ T1 in 0..3,
%@ T2 in 0..3 ;
%@ N1 = N2, N2 = 6,
%@ T1 in 0..3,
%@ T2 in 0..6 ;
%@ T1 = 0,          %% NB: We re-find all same solutions once again...
%@ N1 = N2, N2 = 3,
%@ T2 in 0..3 ;
%@ T1 = 0,
%@ N1 = 3,
%@ N2 = 6,
%@ T2 in 0..6 ;
%@ N1 = 6,
%@ T2 = N2, N2 = 0,
%@ T1 in 0..3 ;
%@ N1 = 6,
%@ N2 = 3,
%@ T1 in 0..3,
%@ T2 in 0..3 ;
%@ N1 = N2, N2 = 6,
%@ T1 in 0..3,
%@ T2 in 0..6.

%% A DESIGN is a RELATION connecting acceptable tallies
%% to subsequent denominators that constitute decisions
%% to enroll more study participants.
oktally_nextdenom(T1/N1 - T2/N2, N1next - N2next) :-
    acceptable_tally(T1/N1 - T2/N2),
    denominator(N1next - N2next),
    N1next #>= N1,
    N2next #>= N2,
    %% At this point, I need to universally quantify the solution
    %% over all possible outcomes.
    E1 #= N1next - N1, % Enrollments at dose 1
    E2 #= N2next - N2, % and dose 2.
    E1 + E2 #> 0, % Exclude trivial enrollment 'steps'
    bagof(T1next/N1next - T2next/N2next,
	  (   tally0_tally(T1/N1 - T2/N2, T1next/N1next - T2next/N2next),
	      indomain(T1next), % What a horror show!
	      indomain(N1next), % There must be a
	      indomain(T2next), % better way!
	      indomain(N2next)
	  ),	      
	  Tallies),
    maplist(acceptable_tally, Tallies).
   
%?- tally0_tally(0/3 - 0/0, T1/6 - 0/0).
%@ T1 in 0..3 ;
%@ T1 in 0..3.

%?- acceptable_tally(0/6 - 0/0).
%@ false.

%?- bagof(T1/6 - 0/0, (tally0_tally(0/3 - 0/0, T1/6 - 0/0), indomain(T1)), Tallies), maplist(acceptable_tally, Tallies).
%@ false.

%?- oktally_nextdenom(0/3 - 0/0, 6 - 0).
%@ false. % CORRECT!

%?- oktally_nextdenom(0/3 - 0/0, 3 - 3).
%@ true ;
%@ true.

%?- oktally_nextdenom(0/3 - 0/0, 3 - N2).
%@ false. %% Non-monotonic behavior! (Blame bagof/2?)

/*
My guess about this is that, if I do want to treat the design as if it were
itself a Prolog program, then I will have to implement some kind of MI.

What can I do, short of this? Perhaps I can simply treat the design as a list
of tuples, in a manner accessible to tuples_in/2. This is perhaps better suited
anyway, to aims such as counting the number of possible designs, etc.

I wonder if I've encountered a specific instance of a general type of tension
that might be common in logic programming, between two conflicting aims. If so,
how could I articulate these aims?
*/

/*
Interesting! What we see here is a set of steps that lead to acceptable tallies.
But what is not checked here is whether the chosen end-denominator ensures that
ALL POSSIBLE end-tallies are acceptable! The latter is a much stronger condition,
and is what in fact yields the USEFUL inference about the CONTENT of a design.

Looking at this from a naive (pre-CLP) perspective, it would seem that I require
the \+ operator. But 2 possibilities present themselves, both anticipated by Markus.

The first is TABLING, which it seems I might use to 'mark' arrows as INVALID after
finding a possible cohort outcome that yields an unacceptable end-tally.

The second is tuples_in/2, which I think would require me to represent tallies
as simple lists of integers.

But TO BEGIN WITH, let me try a very straightforward (even if impure) approach
based on findall/3.
*/


%?- step(T1a/N1a - T2a/N2a, T1b/N1b - T2b/N2b), indomain(N1a), indomain(N2a), indomain(N1b), indomain(N2b).
%@ N1a = N1b, N1b = 6,
%@ T2a = N2a, N2a = T2b, T2b = N2b, N2b = 0,
%@ T1a in 1..4,
%@ T1b in 1..4 ;
%@ N1a = N1b, N1b = 6,
%@ T2a = N2a, N2a = 0,
%@ N2b = 3,
%@ T1a in 1..4,
%@ T1b in 1..4,
%@ T2b in 0..3 ;
%@ N1a = N1b, N1b = 6,
%@ N2a = N2b, N2b = 3,
%@ T1a in 1..4,
%@ T2a in 0..3,
%@ T1b in 1..4,
%@ T2b in 0..3 ;
%@ N1a = N1b, N1b = N2b, N2b = 6,
%@ N2a = 3,
%@ T1a in 1..4,
%@ T2a in 0..3,
%@ T1b in 1..4,
%@ T2b in 0..4 ;
%@ N1a = N2a, N2a = N1b, N1b = N2b, N2b = 6,
%@ T1a in 1..4,
%@ T2a in 0..4,
%@ T1b in 1..4,
%@ T2b in 0..4 ;
%@ N1a = N1b, N1b = 6,
%@ N2a = N2b, N2b = 3,
%@ T1a in 0..4,
%@ T2a in 0..3,
%@ T1b in 1..4,
%@ T2b in 0..3 ;
%@ N1a = N1b, N1b = N2b, N2b = 6,
%@ N2a = 3,
%@ T1a in 0..4,
%@ T2a in 0..3,
%@ T1b in 1..4,
%@ T2b in 0..4 ;
%@ N1a = N2a, N2a = N1b, N1b = N2b, N2b = 6,
%@ T1a in 0..4,
%@ T2a in 0..4,
%@ T1b in 1..4,
%@ T2b in 0..4 ;
%@ N1a = 3,
%@ T2a = N2a, N2a = T2b, T2b = N2b, N2b = 0,
%@ N1b = 6,
%@ T1a in 0..3,
%@ T1b in 1..4 ;
%@ T1a = T2a, T2a = N2a, N2a = T2b, T2b = N2b, N2b = 0,
%@ N1a = 3,
%@ N1b = 6,
%@ T1b in 1..4 ;
%@ T1a = 0,
%@ N1a = N2a, N2a = N2b, N2b = 3,
%@ N1b = 6,
%@ T2a in 0..3,
%@ T1b in 1..4,
%@ T2b in 0..3 ;
%@ T1a = 0,
%@ N1a = 3,
%@ N2a = N1b, N1b = N2b, N2b = 6,
%@ T2a in 0..4,
%@ T1b in 1..4,
%@ T2b in 0..4 ;
%@ N1a = N1b, N1b = 6,
%@ T2a = N2a, N2a = 0,
%@ N2b = 3,
%@ T1a in 1..4,
%@ T1b in 0..4,
%@ T2b in 0..3 ;
%@ N1a = N1b, N1b = 6,
%@ N2a = N2b, N2b = 3,
%@ T1a in 1..4,
%@ T2a in 0..3,
%@ T1b in 0..4,
%@ T2b in 0..3 ;
%@ N1a = N1b, N1b = N2b, N2b = 6,
%@ N2a = 3,
%@ T1a in 1..4,
%@ T2a in 0..3,
%@ T1b in 0..4,
%@ T2b in 0..4 ;
%@ N1a = N2a, N2a = N1b, N1b = N2b, N2b = 6,
%@ T1a in 1..4,
%@ T2a in 0..4,
%@ T1b in 0..4,
%@ T2b in 0..4 ;
%@ N1a = N1b, N1b = 6,
%@ N2a = N2b, N2b = 3,
%@ T1a in 0..4,
%@ T2a in 0..3,
%@ T1b in 0..4,
%@ T2b in 0..3 ;
%@ N1a = N1b, N1b = N2b, N2b = 6,
%@ N2a = 3,
%@ T1a in 0..4,
%@ T2a in 0..3,
%@ T1b in 0..4,
%@ T2b in 0..4 ;
%@ N1a = N2a, N2a = N1b, N1b = N2b, N2b = 6,
%@ T1a in 0..4,
%@ T2a in 0..4,
%@ T1b in 0..4,
%@ T2b in 0..4 ;
%@ T1a = 0,
%@ N1a = N2a, N2a = N2b, N2b = 3,
%@ N1b = 6,
%@ T2a in 0..3,
%@ T1b in 0..4,
%@ T2b in 0..3 ;
%@ T1a = 0,
%@ N1a = 3,
%@ N2a = N1b, N1b = N2b, N2b = 6,
%@ T2a in 0..4,
%@ T1b in 0..4,
%@ T2b in 0..4 ;
%@ T1a = N1a, N1a = T2a, T2a = N2a, N2a = T1b, T1b = N1b, N1b = T2b, T2b = N2b, N2b = 0 ;
%@ T1a = N1a, N1a = T2a, T2a = N2a, N2a = T2b, T2b = N2b, N2b = 0,
%@ N1b = 3,
%@ T1b in 0..3 ;
%@ N1a = N1b, N1b = 3,
%@ T2a = N2a, N2a = T2b, T2b = N2b, N2b = 0,
%@ T1a in 0..3,
%@ T1b in 0..3 ;
%@ T1a = T2a, T2a = N2a, N2a = T2b, T2b = N2b, N2b = 0,
%@ N1a = N1b, N1b = 3,
%@ T1b in 0..3 ;
%@ T1a = N1a, N1a = T2a, T2a = N2a, N2a = T1b, T1b = T2b, T2b = N2b, N2b = 0,
%@ N1b = 3 ;
%@ N1a = N1b, N1b = 3,
%@ T2a = N2a, N2a = T1b, T1b = T2b, T2b = N2b, N2b = 0,
%@ T1a in 0..3 ;
%@ N1a = N1b, N1b = N2b, N2b = 3,
%@ T2a = N2a, N2a = T1b, T1b = 0,
%@ T1a in 0..3,
%@ T2b in 0..3 ;
%@ T1a = T2a, T2a = N2a, N2a = T1b, T1b = T2b, T2b = N2b, N2b = 0,
%@ N1a = N1b, N1b = 3 ;
%@ T1a = T2a, T2a = N2a, N2a = T1b, T1b = 0,
%@ N1a = N1b, N1b = N2b, N2b = 3,
%@ T2b in 0..3 ;
%@ T1a = T1b, T1b = 0,
%@ N1a = N2a, N2a = N1b, N1b = N2b, N2b = 3,
%@ T2a in 0..3,
%@ T2b in 0..3 ;
%@ T1a = T1b, T1b = 0,
%@ N1a = N2a, N2a = N1b, N1b = 3,
%@ N2b = 6,
%@ T2a in 0..3,
%@ T2b in 0..4 ;
%@ T1a = T1b, T1b = 0,
%@ N1a = N1b, N1b = 3,
%@ N2a = N2b, N2b = 6,
%@ T2a in 0..4,
%@ T2b in 0..4.

%% Big questions: what are we talking about at all?

%% We can ABSORB the unification and backtracking!
%% Compare with the call stack of Lisp, which is implicit.

/*

Reading notes from Codish & Søndergaard:

- Various safety properties might lend themselves to analysis in a 'Safety' domain like C&S's 'Parity' domain.

- The abstract elements 'even', 'odd', etc., look like 'escalation', 'de-escalation', and other such 'commentary' you might offer observing the enrollment and assessment of cohorts.

- Abstract interpretation as described in introduction to Section 3 looks like a method for demonstrating higher-level concepts as 'emergent' from the designs as demarcated by merely 'algorithmic' constructs.

- TODO: Try to understand if/how 'fixed points' relate to the above construction of 'designs'.

- In practical terms, how would a Tp semantics for dose-escalation trials differ from a 'standard meta-circular interpreter'? What does the symmetric difference of the sets of capabilities of these approaches look like? (I'm asking, what I may appreciate only in retrospect, how the considerations of Section 3.1 manifest concretely in this particular application.)

- Is 'least fixed point' something like a *transitive closure* in RDBMS?

- What's missing from the fixed point of append/3 as given on page 8? The possibility of 'creative' uses such as with list differences?

- DISCUSS: If one does pure logic programming, does this eliminate (except say for purely technical reasons like debugging or performance) all interest in "runtime" issues such as in Section 3.3?

- I could see depth-k analysis achieving a 'local' form of dose-escalation trial analysis that focuses attention on (say) the last move, present decision, and next set of outcomes/decisions. Some of the discussion initiated by Wages & Braun has this this character.

- Does the "irrational assignment" question of Wages & Braun link to 'dataflow analysis'? You're asking the question, whether the outcome for the next patient changes ('flows through to') the future dosing decisions. 

- I even wonder whether the DCG {esc,sta,des} amounts to an 'approximation' or could be recast as such (compare Section 4.3, where predicates are approximated by terms).
A: But esc/sta/des are already the finest-grained things that can happen. Contrast with even/odd. When you say "approximation" you're talking about domain granularity.

- Approximation within the 3+3 design could make explicit the equal treatment given to 2/3 and 3/3 cohorts -- and even pave the way for a *criticism* of this equivalence! (A similar point might well apply to representing and criticizing the dichotomization of ordinal (Gr in 0..5) toxicities via DLT in 0..1).
A: But how much can we abstract from details, yet retain useful consequences that we can 'milk' from it.

- Would it help me to know the origins of 's' and 'c' in (s|c)-semantics?
TODO: Look up the papers (e.g. ref. 7, 30, 31).

- The manner in which MI's abstract/extend search strategies seems potentially like a heuristic or metaphor for what one would like to achieve wrt dose-escalation designs. There are in general some 'absorbed' notions that remain unchanged, but other 'reified' notions that admit modification.

*/

/*
TODO: Consider metagol!
Homoiconicity.
Functor-free subset of Prolog is Datalog!
*/
