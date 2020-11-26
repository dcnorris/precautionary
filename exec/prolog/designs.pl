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
    maplist(indomain, [E1, E2]),
    %indomain(E1), % Avoid 'arguments not sufficiently
    %indomain(E2), % instantiated' error
    T1new in 0..E1, % TODO: Try plainer constraints here,  ***
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
    %% TODO: Why can't MAPLIST do the same thing BY PURE MEANS?
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
%@ false. %% NON-MONOTONIC BEHAVIOR! (Blame \+ via bagof/2?)

/*
My guess about this is that, if I do want to treat the design as if it were
itself a Prolog program, then I will have to implement some kind of MI.

What can I do, short of this? Perhaps I can simply treat the design as a list
of tuples, in a manner accessible to tuples_in/2. This is perhaps better suited
anyway, to aims such as counting the number of possible designs, etc.

I wonder if I've encountered a specific instance of a general type of tension
that might be common in logic programming, between two conflicting aims. If so,
how could I articulate these aims? Is it a distinction between REIFICATION and
a more 'subtle' or 'invisible' or IMPLICIT representation?

Or perhaps the question is merely that of TEMPORAL vs SPATIAL REPRESENTATION?
*/

%% I feel as if I'm in need of GENERAL advice about 'what to do when you feel
%% like using an all-solutions predicate'. The alternatives that come to mind
%% include \+ and cut -- which at least are more honest about how bad they are.

/*
Maybe this situation forces me to employ CLP in a more essential way.

As unfortunate as tuples_in/2's (unstructured) lists of integers might appear,
they do promise a certain ease of generalizability toward D > 2.
The transformation from tallies like T1/N1 - T2/N2 to tuples [T1, N1, T2, N2]
isn't really so very hard to read, either!
*/


%% Suppose I tried to express the constraint via tuples_in/2?
/*
What about tmember/2, tfilter/3 and tpartition/4?
*** These are in library(reif).
*/

%% Big questions: what are we talking about at all?

%% We can ABSORB the unification and backtracking!
%% Compare with the call stack of Lisp, which is implicit.

/*

Residual notes from C&S reading:

- I could see depth-k analysis achieving a 'local' form of dose-escalation trial analysis that focuses attention on (say) the last move, present decision, and next set of outcomes/decisions. Some of the discussion initiated by Wages & Braun has this this character.

- Does the "irrational assignment" question of Wages & Braun link to 'dataflow analysis'? You're asking the question, whether the outcome for the next patient changes ('flows through to') the future dosing decisions. 

- The manner in which MI's abstract/extend search strategies seems potentially like a heuristic or metaphor for what one would like to achieve wrt dose-escalation designs. There are in general some 'absorbed' notions that remain unchanged, but other 'reified' notions that admit modification.

*/
%% Decouple the description frrom the DCG formalism.
%% Free the description.
%% TODO: Inspect Markus's 3 DSLs from clpz.
%% "Finding out is always an example of search." -- MT

/*
Notes from reading library(clpz) ...

WHAT IF there turned out to be some clever encoding of 'designs'
in terms of integer and combinatorial constraints?

Perhaps some kind of (conditional?) partial ordering is conceivable,
which would avoid the \+ problem? Another way to pose this question is:
what conditions would have to hold for our 'unacceptabilty' notions,
in order to ensure that tallies may be mapped to Z to render the
acceptable sets as pre-images of INTERVALS in Z?

-) I think the tallies over a given denominator might well be amenable
   to a complete ordering -- on the basis of likelihood of toxicity.
   If this proved to be the case, then this would constitute a marvelous
   linkage to the analytical or topological content of pharmacology in
   the continuum ℝ.
-) If this mapping can be done bijectively, then I suddenly have a problem
   entirely within the capabilities of CLP(ℤ).
-) Even if no bijection of this kind can be stated a priori, perhaps we may
   determine one 'dynamically', consistently with whatever acceptability
   constraints have been imposed by the domain expert.
-) What I may require, in order to render this possible, is some arithmetical
   fact---with luck, a TRUE one!---about tallies in general. Do they admit an
   ordering that renders all possible outcomes of any given cohort enrollment
   as an interval within the order? Does this ordering enable any acceptability
   criterion to be expressed as an interval?
*/

/*
ADDITIONAL THOUGHTS ..

What about that recursive (dynamic programming) perspective I had in mind a while back?
Would backwards connecting eliminate the need for forwards universal quantification?
*/

/* - - - - -

What if I revisit my LATTICE concept, with the assistance of an
encoding that brings CLP(Z) more effectively into play?

Let's also be less shy about using LOGICAL VARIABLES (FTW!),
so that our DENOMINATOR LATTICE looks like this:

6  I..G..E
   .......
   .......
3  H..F..D
   .......
   .......
0  A..B..C
   0  3  6

In part, the above ordering corresponds to lexicographic preference
for escalating 'as little and as late' as possible. But I have also
aimed to place the final letters {G, H, I} at the denominators least
likely (or impossible) to reach by cautious dose-escalation steps.

Can I formalize this intuition? Writing the denominators low-dose-last,
I obtain this layout:

6  60..63..66
   ..........
   ..........
3  30..33..36
   ..........
   ..........
0  00..03..06
    0   3   6

Interestingly, while I can see no purely arithmetical treatment of
the denominators *themselves* which represents my intuition, I do
find that PATH concepts help with the ordering. The letters A..E
trace out a path of maximally-delayed dose escalation. The next path
in that ordering is ABFDE, and the next after that ABFGE. Thus, the
labeling of the nodes is such that, when listed in order of ascending
'eagerness', the escalation paths introduce the later letters of the
alphabet at the latest time possible. (The first path to 66 employs
only the first 5 alphabet letters; the next path, the first 6 letters,
and so on.)

Notice also that the integer encoding of the denominators DOES nicely
indicate the 'eagerness' of the paths, if the integers themselves are
treated lexicographically:

(00,03,06,36,66)
(00,03,33,36,66)
(00,03,33,63,66)
(00,30,33,36,66)
(00,30,33,63,66)
(00,30,60,63,66).

Those are the 6 possible dose-escalation paths running from 00 to 66
that are STRICTLY INCREASING in the XY encoding. (Note that 6 is the
binomial coefficient (4 choose 2), and that the path-counting here
recapitulates Pascal's triangle.)

Now that I'm not afraid to create 9 logical variables from the outset,
what becomes possible in the way of posting constraints?

First of all, can I exploit 'meta-logical' predicates, such that I
ask whether VARIABLES THEMSELVES can be linked by escalation steps?
If I can 'get away with' this, then maybe I have an escape from \+?

Here's the proposed set-up: The variables A..I are bound to lists of
numerators, coded in the XY encoding as above. Two variables V1, V2
may be linked by a valid dose-escalation step iff there is a numerator
in the list bound to V1 that has a NONDETERMINISTIC IMAGE entirely
within the list bound to V2.

Notice that this is a WEAKER condition than that imposed on the
'escalation arrows' relation conceived above. The variables are
connected if there is ANY SUCH ARROW. Thus, the 'links' I am
discussing here are equivalence classes of a sort.

===

Getting into SPECIFICS now, here's how these arrows might be tested.

Say B = [00,01,02,03] to begin with. As it turns out, none of these
is (or, logically, even *could* be) unacceptable. We now ask what
arrows might be drawn:

00 -> 06 (C) x should escalate
00 -> 33 (F)
01 -> 06 (C)
01 -> 33 (F) x should stay
02 -> 06 (C) x should stop
02 -> 33 (F) x should stop
03 -> 06 (C) x should stop
03 -> 33 (F) x should stop

This is the full set of arrows from B-POSSIBLE TALLIES to one of
the feasible target denominators C ('06') or F ('33'). Disallowed
arrows are marked 'x', with an explanation.

Can I exhibit CLP(Z)-type constraints that yield the x's above?
Generating the DISALLOWED arrows ought to work just fine, provided
that I am then willing to move to a SPATIAL representation over FD,
and perform the necessary set-differencing.

OOH! But there is clearly a problem with the feasibility of generating
enormous solution sets, only to difference them to obtain small ones!
Really, the first question I ought to ask is whether I can obtain the
ALLOWED arrows via CLP(Z) constraint propagation.

Perhaps I should focus for a while on the problem of achieving a
representation of the INDETERMINACY itself. Certainly, CLP(Z) offers
handy forms of expression, via .. , /\ and \/ notations.

Do any CLP(Z) combinatorial constraints come at all close to 'forall'?

-) It seems at least that all_distinct/1 does an awful lot of work!
-) What about disjoint2(+Rectangles) ...

AHA! This RECTANGLES idea makes me wonder if I have neglected to seek
out SPATIAL INTUITIONS for this problem! Do the sets of possible cohort
outcome tallies have geometrical analogues? Between any two denominators
we can draw a rectangle, and the possible cohort tallies lie in the
origin-translated image of that rectangle. What about the end-tallies
themselves? These numerators obey the same rules. Can my rules about
(un)acceptablity be translated to the origin as well? Are convexities
of any sort thereby revealed? 

Anticipating (for a moment) the possibilities contained in disjoint2/1,
once the UNacceptable rectangles have been mapped out in the tally grid,
I ought to be able to identify acceptable escalations as rectangles that
are disjoint from the union of all the unacceptable regions. (Then again,
maybe this will look like overkill, with the advantage of sufficient
geometrical intuition?)

===
What are the REGRETS in the table above?

00 -> 06 (C) x 06/06
00 -> 33 (F)
01 -> 06 (C)
01 -> 33 (F) x 31/33, 21/33
02 -> 06 (C) x 05/06
02 -> 33 (F) x 32/33 (and others?)
03 -> 06 (C) x 05/06
03 -> 33 (F) x Hmm...

Ah, I see that not every unacceptable tally is defined purely by a
regret about the net end result. Maybe this is a fault in my present
formulation! Perhaps REJECTING THE DECISIONS THEMSELVES 'gets ahead
of myself' a bit, and the more coherent approach will be based upon
true REGRETS ABOUT OUTCOMES. (Under this latter approach, criticism
of the decisions in themselves is properly seen as a 2nd-order idea;
the primary ideas should perhaps be our regrets about net tallies,
irrespective of how we reached them.)

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



   - - - - - */
?- clpz:must_be(X, Y, a).
%@    X = ground.

?- clpz:must_be(acyclic, Y, a).
%@    true.
%@ caught: error(type_error(type,acylic),must_be/2)
%@ caught: error(type_error(type,acylic),must_be/2)
%@ caught: error(evaluation_error((clpz:must_be)/3),must_be/3)

:- use_module(library(clpz)).
%@    true.
