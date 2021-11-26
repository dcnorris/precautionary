%% At LONG last I feel ready to sketch a DSL in this domain ...

/*

Here is a path-probability-weighted summary of des/sta/esc decisions
happening during operation of the Braun (2020) CRM model with C=13.
Note that there is a pretty sharp demarcation along boundaries in the
T-N plane.

This inspires for me the hope that designs formulated via constraints
may closely approximate even CRM designs.

, , E = des

   N
T   1         2 3           4 5           6 7            8 9           10
  0 0 0.0000000 0 0.000000000 0 0.000000000 0 0.0000000000 0 0.0000000000
  1 0 0.1822642 0 0.023113476 0 0.000000000 0 0.0000000000 0 0.0000000000
  2 0 0.2118588 0 0.181816722 0 0.048386476 0 0.0010110676 0 0.0001071718
  3 0 0.0000000 0 0.094821725 0 0.146228250 0 0.0304239336 0 0.0060178019
  4 0 0.0000000 0 0.006561661 0 0.031121860 0 0.0680561438 0 0.0696644758
  5 0 0.0000000 0 0.000000000 0 0.003347426 0 0.0143452867 0 0.0302239098
  6 0 0.0000000 0 0.000000000 0 0.000000000 0 0.0006738165 0 0.0043339708

, , E = sta

   N
T   1         2 3          4 5            6 7            8 9          10
  0 0 0.1841717 0 0.28902885 0 8.707279e-02 0 2.240090e-02 0 0.014560425
  1 0 0.7317945 0 0.66745341 0 3.092600e-01 0 2.309519e-01 0 0.144591290
  2 0 0.0009000 0 0.19450867 0 3.622027e-01 0 3.697338e-01 0 0.298500644
  3 0 0.0000000 0 0.00010476 0 2.659965e-02 0 2.016770e-01 0 0.256227656
  4 0 0.0000000 0 0.00000081 0 6.530033e-03 0 1.531475e-02 0 0.068165189
  5 0 0.0000000 0 0.00000000 0 1.414260e-07 0 1.074698e-06 0 0.009602811
  6 0 0.0000000 0 0.00000000 0 7.290000e-10 0 0.000000e+00 0 0.000000000

, , E = esc

   N
T   1        2 3          4 5         6 7          8 9          10
  0 0 2.246438 0 0.18293632 0 0.1606449 0 0.07546003 0 0.020376738
  1 0 0.000000 0 0.05476038 0 0.2162357 0 0.08066563 0 0.034548882
  2 0 0.000000 0 0.00000000 0 0.0000000 0 0.00000000 0 0.002925861
  3 0 0.000000 0 0.00000000 0 0.0000000 0 0.00000000 0 0.000000000
  4 0 0.000000 0 0.00000000 0 0.0000000 0 0.00000000 0 0.000000000
  5 0 0.000000 0 0.00000000 0 0.0000000 0 0.00000000 0 0.000000000
  6 0 0.000000 0 0.00000000 0 0.0000000 0 0.00000000 0 0.000000000

> 
*/


/*

I would like to specify dose-finding intuitions in the form of
forbidden events. With a suitable interpreter, I believe that
the 3+3 trial will lend itself to such definition.

Leaving aside for the moment the question of exactly HOW the
interpreter should work, let me articulate some constraints
that /with suitable interpretation/ MAY yield the 3+3 design.

===

Q1: Why do we escalate after observing 0/3?

A1: Because regret([sta, 0/6]).

---

Q2: Why do we DEescalate after 2+/3?

A2: Because regret([sta, 5/6]).

---

Q3: Why do we stay after 1/3?

A3: Hmm.. This one might be the source of some insight!

- What if I escalated after 1/3, and saw 3/3? Is this
  the outcome that looks regrettable in retrospect?
  How do I /articulate/ this regrettableness?

- OR.. is this better expressed in terms other than 'regret'?
  That is, should I say that escalating requires some 'support',
  which an observation of 0/3 provides?

- Are there HIGHER-LEVEL concepts at work, that deserve
  to be expressed in the DSL? For example, are there cohort
  sizes that I say inadequately explore a dose FOR A GIVEN
  PURPOSE?

- Then again, the excessive proliferation of extra concepts
  should be avoided assiduously! In these early explorations,
  my very first question should be, How much can I accomplish
  with a very small language?

- One useful principle might be that we always choose the next
  dose to be as high as allowed. (Maybe this is even a parameter
  that we could change via DSL!) If so, then Q3 gets further
  focused as, "Why DON'T we escalate after 1/3?"

- Keep in mind that Q3 may provoke the articulation of SEVERAL
  distinct principles. Remember, my aim is to develop insight
  at this stage!

===

Am I perhaps aiming at a PROGRAMMING-BY-EXAMPLE principle?
No, it's https://en.wikipedia.org/wiki/Inductive_programming!

*/

%% I can start the programming, already!
:- use_module(library(lists)).
:- use_module(library(clpz)).
:- use_module(library(reif)).
:- use_module(library(si)).
:- use_module(library(dcgs)).
:- use_module(library(error)).
:- use_module(tally).

%% Switching to a tidier representation of state that does not
%% carry along excluded doses...
%% Also, adopting the convention that the current dose is the
%% head of the left-hand list.

enroll(T0/N0, T1/N1) :-
    #Nnew #= 3, % hard-coding cohorts of 3 for now
    #N1 #= #N0 + #Nnew,
    Tnew in 0..Nnew, indomain(Tnew),
    #T1 #= #T0 + #Tnew,
    #N1 #=< 6. % Maximum enrollment per cohort
    
%?- enroll(0/0, T/n).
%@ caught: error(existence_error(procedure,enroll/2),enroll/2)
%@ caught: error(existence_error(procedure,enroll/2),enroll/2)

state0_decision_state(Ls - [R0|Rs], esc, [R|Ls] - Rs) :-
    enroll(R0, R).

state0_decision_state([L0|Ls] - Rs, sta, [L|Ls] - Rs) :-
    enroll(L0, L).

state0_decision_state([L,D0|Ls] - Rs, des, [D|Ls] - [L|Rs]) :-
    enroll(D0, D).

%?- state0_decision_state([0/3, 1/6] - [0/0, 0/0], esc, Ls - Rs).
%@    Ls = [0/3,0/3,1/6], Rs = [0/0]
%@ ;  Ls = [1/3,0/3,1/6], Rs = [0/0]
%@ ;  Ls = [2/3,0/3,1/6], Rs = [0/0]
%@ ;  Ls = [3/3,0/3,1/6], Rs = [0/0]
%@ ;  false.

%?- state0_decision_state([1/3, 1/6] - [0/0, 0/0], sta, Ls - Rs).
%@    Ls = [1/6,1/6], Rs = [0/0,0/0]
%@ ;  Ls = [2/6,1/6], Rs = [0/0,0/0]
%@ ;  Ls = [3/6,1/6], Rs = [0/0,0/0]
%@ ;  Ls = [4/6,1/6], Rs = [0/0,0/0]
%@ ;  false.

%?- state0_decision_state([2/3, 0/3] - [0/0, 0/0], des, Ls - Rs).
%@    Ls = [0/6], Rs = [2/3,0/0,0/0]
%@ ;  Ls = [1/6], Rs = [2/3,0/0,0/0]
%@ ;  Ls = [2/6], Rs = [2/3,0/0,0/0]
%@ ;  Ls = [3/6], Rs = [2/3,0/0,0/0].

%% What DCG could describe valid sequences of (Action, Tally) pairs?
%% Is a list of such pairs really the right representation? Might a
%% fluid, alternating sequence of tallies & decisions work better?

%?- state0_decision_state([] - [0/0, 0/0], esc, S).
%@    S = [0/3]-[0/0]
%@ ;  S = [1/3]-[0/0]
%@ ;  S = [2/3]-[0/0]
%@ ;  S = [3/3]-[0/0]
%@ ;  false.

/*

PLAN from morning walk 8/28 ...

1. Generate common designs from regret-constrained (RC) formulations
(a) Standard 3+3
(b) Less common variant
(c) CCDs   

2. Rolling enrollment
(a) 3+3 à la Frankel &al (2020) /JAMA Open/
(b) Rolling-6 (Skolnik &al 2008)

3. Approximating CRM
(a) Exhibit an RC design that contains all paths of the Braun2020 model.
    Demonstrate how much larger the path set is (probably a lot!) and
    characterize what decisional indeterminacy remains in the design.
(b) Find a principled way to remove sets of paths by adding RCs, and
    to characterize the path-set differences that arise along the way
    toward a design that is contained within the CRM.

*/

% TODO:
% 1. How do I exclude possibilities without not/1 or all-solutions predicates?
%
% What is the proper STATUS of statements like regret/2?
% Are these actually to be implemented differently?

/*

I need to think 'from a logical point of view' at this point.
I would like to posit EXAMPLES of regrettable situations, and have these RESTRICT
the possible solutions. But this runs counter to the way pure logical programming
is supposed to work, such that adding more facts should only increase the number
of solutions.

Thus, the statement of regrettable situations is more like specifying CONSTRAINTS
as opposed to FACTS.

This reasoning seems to demand that I specify 'regret' in terms of CLP, and not
in terms of predicates or facts. The design is not written in Prolog, but in CLP!

*/

/*

I am beginning to think that the above observation is the crucial one. If I keep
trying to handle regret using Prolog, I'm working at cross purposes.

But now the question becomes, how do I incorporate CLP statements of regret into
the operation of the DCG that generates possible trial paths?

----------------------------------------

WAIT! Why can't I add 'regret' statements to Prolog's database, as constraints?
Why can't Prolog search every possible reason for regretting an action, and then
if it can't find one, do the action? Does this REALLY require not/1?

I THINK IT DOES! The trouble is, unless I commit to the discovery of a regret,
upon backtracking the list of regrets will get exhausted until the regretted
action is eventually taken anyway! Thus, taking the action requires satisfying
the goal of *NOT* REGRETTING. And there's the not/1!

BUT... might if_/3 be able to rescue me? Does it discard choice points purely?
If I understand it, YES!

But the price to pay for this is REIFICATION. I need to *reify* regret.

*/

%% regret(A, [Q|Qs]) holds if we regret having taken action A if it
%% results in a tally history [Q|Qs]. Note that this allows, in practice,
%% for the relevant history to go back 1 or 2 tallies.
%% TODO: Consider whether this list 'trick' is a bad thing. Should I just
%%       create an extra argument, invoking with _ in don't-care cases?

regret(sta, [T/N|_]) :-
    (	#N #> 6; % regret enrolling too many at any one dose
	#N #= 6 #/\ #T #= 0; % regret staying at a too-low dose
	#N #= 6 #/\ #T #>= 5 % regret excess toxitiy
    ).

%?- regret(sta, [T/N|_]).
%@    clpz:(N in 7..sup)
%@ ;  T = 0, N = 6
%@ ;  N = 6, clpz:(T in 5..sup).

%?- regret(sta, [5/6]).
%@    true.

regret(esc, [T/3, T0/3]) :- T in 0..3, T0 in 0..3,
			    #T0 #> 0 #/\ #T #= 3.

%% This clause expresses that we regret ANY toxicities after
%% having escalated from a dose with fewer than N=3 assessments.
regret(esc, [T/3, T0/N0]) :- T in 0..3, N0 in 0..6, T0 in 0..6,
			     #T0 #=< #N0,
			     #T #> 0 #/\ #N0 #< 3.

%% The above is not quite sufficient, however, since we also
%% ought to regret ANY toxicities after escalating from >= 2/6.
regret(esc, [T/N, T0/6]) :- T0 in 0..6, N in 0..6, (#N #= 3 #\/ #N #= 6), T in 0..6, #T #=< #N,
			    #T #> 0,
			    #T0 #> 1.

%?- regret(esc, [T/N, T0/N0]).
%@    T = 3, N = 3, N0 = 3, clpz:(T0 in 1..3)
%@ ;  N = 3, clpz:(#N0+1#= #_A), clpz:(#N0#>= #T0), clpz:(T in 1..3), clpz:(_A in 1..3), clpz:(N0 in 0..2), clpz:(T0 in 0..2)
%@ ;  N0 = 6, clpz:(#N#=3#<==> #_A), clpz:(#N#=6#<==> #_B), clpz:(#_A#\/ #_B#<==>1), clpz:(#N#>= #T), clpz:(_A in 0..1), clpz:(_B in 0..1), clpz:(N in 1..6), clpz:(T in 1..6), clpz:(T0 in 2..6).

%?- regret(esc, Qs).
%@    Qs = [3/3,_A/3], clpz:(_A in 1..3)
%@ ;  Qs = [_D/3,_C/_A], clpz:(#_A+1#= #_B), clpz:(_C in 0..sup), clpz:(_D in 1..3), clpz:(_B in 1..3), clpz:(_A in 0..2)
%@ ;  Qs = [_B,_A/6], clpz:(_A in 2..6).

%% This holds when decision A is a dose-escalation decision
%% that is FEASIBLE from state S0, and for which 'regret'
%% is an applicable concept. The Truth-value reifies whether
%% the action might (under some realization of toxicities)
%% be judged regrettable, or not. Note that this goal simply
%% fails when A is infeasible ('inconceivable') from S0,
%% for example when A=esc and S0=_-[] (such that there is no
%% higher dose to escalate into).
state0_decision_regrettable(S0, A, Truth) :-
    state_si(S0),
    member(A, [esc,sta]), % these are the decisions to which 'regret' applies
    (	state0_decision_state(S0, A, S),
	S0 = [T0/N0|_] - _, % TODO: Factor this pattern matching
	S  = [T /N |_] - _, %       into a regret/3 predicate?
	%% I introduce the (->) below to avert backtracking over
	%% possibly multiple scenarios for regret -- one is enough!
    	regret(A, [T/N, T0/N0]) -> Truth = true
    ;	state0_decision_state(S0, A, S) -> % ensure A is feasible, but don't re-check this
	Truth = false
    ).

%% I believe the semantics of 'regret' are untenable at the level of if_/3
%% application. What I need is a goal that reifies to 'true' for precisely
%% those decisions which are BOTH feasible AND non-regrettable, and otherwise
%% reifies to 'false'.
state0_decision_noregrets(S0, A, Truth) :-
    state_si(S0),
    member(A, [esc,sta]), % these are the decisions to which 'regret' applies
    %% TODO: Are there advantages of if_/3 even at this level?
    %%       Or would I just start an infinite regress? Is there
    %%       always a '->' at some level (as in if_/3 definition!)
    %%       in any code doing this sort of thing?
    (	state0_decision_state(S0, A, _) -> % commit if feasible
	(   state0_decision_state(S0, A, S),
	    S0 = [T0/N0|_] - _, % TODO: Factor this pattern matching
	    S  = [T /N |_] - _, %       into a regret/3 predicate?
	    %% I introduce the (->) below to avert backtracking over
	    %% possibly multiple scenarios for regret -- one is enough!
    	    regret(A, [T/N, T0/N0]) -> Truth = false
	;   Truth = true
	)
    ;	Truth = false % from S0, A is not feasible
    ).

%?- state0_decision_noregrets([0/3,1/6]-[], E, Truth).
%@    E = esc, Truth = false
%@ ;  E = sta, Truth = false
%@ ;  false.

%% TODO: Does this (if-then) count as idiomatic usage of list_si/1?
%%       Without it, I'm left with unsightly extra choice points.
state_si(L - R) :-
    %% TODO: Should I also insist that all entries in L & R
    %%       are tallies T/N, with T and N being integers?
    %%       (NB: list_si/1 does sort/2, which implies visiting
    %%       every element of the list anyway.)
    (	list_si(L) -> list_si(R)
    ;	false
    ).

%?- state0_decision_regrettable([1/3]-[0/0], E, T).
%@    E = esc, T = true
%@ ;  E = sta, T = false
%@ ;  false.

%?- state0_decision_noregrets([1/3]-[0/0], E, T).
%@    E = esc, T = false
%@ ;  E = sta, T = true
%@ ;  false.

%?- state0_decision_regrettable([0/3]-[0/0], E, T).
%@    E = esc, T = false
%@ ;  E = sta, T = true
%@ ;  false.

%?- state0_decision_regrettable([1/3]-[0/0, 0/0], esc, Truth).
%@    Truth = true.

%% Even if the nesting here feels a bit difficult to read, this code reflects
%% the primacy of REGRET as the crucial user-level concept shaping these designs.

%path(_) --> []. % a convenience for testing; path can stop at any time
path(declare_mtd(_)) --> [].
path(S0) --> { if_(state0_decision_noregrets(S0, esc),
		   (   E = esc,
		       state0_decision_state(S0, E, S)
		   ),
		   if_(state0_decision_noregrets(S0, sta),
		       (   E = sta,
			   state0_decision_state(S0, E, S)
		       ),
		       (   state0_decision_state(S0, des, _) ->
			   (   E = des,
			       state0_decision_state(S0, E, S)
			   )
		       ;   E = stop, % ..and otherwise, STOP.
			   MTD = todo, % TODO: Actually obtain MTD as integer >= 0.
			   S = declare_mtd(MTD)
		       )
		      )
		  )
	     },
	     [E, S],
	     path(S). % TODO: Implement declare_mtd possibility

%?- Path = [sta,[3/3]-[0/0],stop,declare_mtd(todo)], phrase(path([0/0]-[0/0]), Path).
%@    Path = [sta,[3/3]-[0/0],stop,declare_mtd(todo)]
%@ ;  false.

%?- phrase(path([]-[0/0,0/0]), Path).
%@    Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[0/6]-[0/3],esc,[0/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[0/6]-[0/3],esc,[1/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[0/6]-[0/3],esc,[2/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[0/6]-[0/3],esc,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[1/6]-[0/3],esc,[0/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[1/6]-[0/3],esc,[1/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[1/6]-[0/3],esc,[2/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[1/6]-[0/3],esc,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[2/6]-[0/3],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[],des,[3/6]-[0/3],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[1/6,0/3]-[],des,[0/6]-[1/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[1/6,0/3]-[],des,[1/6]-[1/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[1/6,0/3]-[],des,[2/6]-[1/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[1/6,0/3]-[],des,[3/6]-[1/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[0/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[1/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[2/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[3/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[0/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[1/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[2/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[3/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[0/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[1/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[2/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[3/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[0/6]-[2/3],esc,[2/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[0/6]-[2/3],esc,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[0/6]-[2/3],esc,[4/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[0/6]-[2/3],esc,[5/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[1/6]-[2/3],esc,[2/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[1/6]-[2/3],esc,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[1/6]-[2/3],esc,[4/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[1/6]-[2/3],esc,[5/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[2/6]-[2/3],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[3/6]-[2/3],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[0/6]-[3/3],esc,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[0/6]-[3/3],esc,[4/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[0/6]-[3/3],esc,[5/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[0/6]-[3/3],esc,[6/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[1/6]-[3/3],esc,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[1/6]-[3/3],esc,[4/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[1/6]-[3/3],esc,[5/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[1/6]-[3/3],esc,[6/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[2/6]-[3/3],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[3/6]-[3/3],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[1/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[2/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[4/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[2/3,1/6]-[],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[3/3,1/6]-[],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[2/6]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[3/6]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[1/3]-[0/0],sta,[4/6]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[2/3]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [esc,[3/3]-[0/0],stop,declare_mtd(todo)]
%@ ;  false.

%?- state0_decision_state([1/3]-[0/0], sta, S).
%@    S = [1/6]-[0/0] % the T=0 possibility does get 'seen'
%@ ;  S = [2/6]-[0/0]
%@ ;  S = [3/6]-[0/0]
%@ ;  S = [4/6]-[0/0]
%@ ;  false.

%?- state0_decision_noregrets([1/3]-[0/0], E, T).
%@    E = esc, T = false
%@ ;  E = sta, T = true
%@ ;  false.

%?- state0_decision_regrettable([0/3,1/6]-[], E, T).
%@    E = sta, T = true
%@ ;  false.

%% AHA! HERE is the problem. The conditional clause of if_/3 must
%% actually succeed (reifying to 'true' or 'false') in order for
%% the if_/3 goal to succeed.
%% By creating a 3-valued logic here, I have violated this design.
%% So I need a different semantics!
%?- state0_decision_regrettable([0/3,1/6]-[], esc, T).
%@ false.

%?- state0_decision_state([0/3,1/6]-[], E, S).
%@    E = sta, S = [0/6,1/6]-[]
%@ ;  E = sta, S = [1/6,1/6]-[]
%@ ;  E = sta, S = [2/6,1/6]-[]
%@ ;  E = sta, S = [3/6,1/6]-[]
%@ ;  E = des, S = [1/9]-[0/3]
%@ ;  E = des, S = [2/9]-[0/3]
%@ ;  E = des, S = [3/9]-[0/3]
%@ ;  E = des, S = [4/9]-[0/3].

%% AHA. Once again, I'm caught calling a decision regrettable
%% even when it is not feasible.

%?- state0_decision_regrettable([1/6]-[0/0], E, T).
%@    E = esc, T = false
%@ ;  E = sta, T = true
%@ ;  false.

%?- state0_decision_state([1/6]-[0/0], esc, S).
%@    S = [0/3,1/6]-[]
%@ ;  S = [1/3,1/6]-[]
%@ ;  S = [2/3,1/6]-[]
%@ ;  S = [3/3,1/6]-[]
%@ ;  false.

%% So the problem is a failure to stop appropriately from some states?

%?- phrase(path([2/3,1/6]-[]), Path).
%@    Path = []
%@ ;  false.

%?- state0_decision_regrettable([2/3,1/6]-[], esc, true).
%@ false.

%?- state0_decision_state([2/3,1/6]-[0/0], E, S).
%@    E = esc, S = [0/3,2/3,1/6]-[]
%@ ;  E = esc, S = [1/3,2/3,1/6]-[]
%@ ;  E = esc, S = [2/3,2/3,1/6]-[]
%@ ;  E = esc, S = [3/3,2/3,1/6]-[]
%@ ;  E = sta, S = [2/6,1/6]-[0/0]
%@ ;  E = sta, S = [3/6,1/6]-[0/0]
%@ ;  E = sta, S = [4/6,1/6]-[0/0]
%@ ;  E = sta, S = [5/6,1/6]-[0/0]
%@ ;  E = des, S = [1/9]-[2/3,0/0]
%@ ;  E = des, S = [2/9]-[2/3,0/0]
%@ ;  E = des, S = [3/9]-[2/3,0/0]
%@ ;  E = des, S = [4/9]-[2/3,0/0].

%?- state0_decision_state([1/6]-[0/0], E, S).
%@    E = esc, S = [0/3,1/6]-[]
%@ ;  E = esc, S = [1/3,1/6]-[]
%@ ;  E = esc, S = [2/3,1/6]-[]
%@ ;  E = esc, S = [3/3,1/6]-[]
%@ ;  E = sta, S = [1/9]-[0/0]
%@ ;  E = sta, S = [2/9]-[0/0]
%@ ;  E = sta, S = [3/9]-[0/0]
%@ ;  E = sta, S = [4/9]-[0/0]
%@ ;  false.

%?- state0_decision_regrettable([1/6]-[0/0], E, true).
%@    E = sta.

%?- state0_decision_regrettable([1/6]-[0/0], esc, Truth).
%@    Truth = false.

%% TODO: Repair this MGQ!
%?- state0_decision_regrettable([1/6]-[0/0], E, Truth).
%@    E = sta, Truth = true
%@ ;  false. % Lost solution (E=esc, Truth=false)!

%% Aha! So I haven't coded regret for 'des' correctly.
%?- state0_decision_regrettable([3/3]-[0/0], des, Truth).
%@    Truth = false.

%% Or maybe the problem is that my initial intuition about the non-regrettability of 'des'
%% was on the right track somehow. Whereas for 'esc' and 'sta' regrettability requires that
%% the action be possible, this isn't quite what I need at the very bottom of this if_/3
%% cascade.
%% As the last dose-escalation decision possible, 'des' has to be 'regrettable' for slightly
%% different reasons -- which may not even be properly termed 'regret'! Here, we need a 'true'
%% reification where 'des' is EITHER impossible OR 'regrettable'.

%?- length(Path, L), phrase(path([0/0]-[0/0,0/0]), Path).
%@ 
%@ caught: error('$interrupt_thrown',repl)
%@    Path = [], L = 0
%@ ;  Path = [sta,[0/3]-[0/0,0/0]], L = 2
%@ ;  Path = [sta,[1/3]-[0/0,0/0]], L = 2
%@ ;  Path = [sta,[2/3]-[0/0,0/0]], L = 2
%@ ;  Path = [sta,[3/3]-[0/0,0/0]], L = 2
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[0/3,0/3]-[0/0]], L = 4
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0]], L = 4
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0]], L = 4
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0]], L = 4
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[1/6]-[0/0,0/0]], L = 4
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[2/6]-[0/0,0/0]], L = 4
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[3/6]-[0/0,0/0]], L = 4
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[4/6]-[0/0,0/0]], L = 4
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[0/3,0/3]-[0/0],esc,[0/3,0/3,0/3]-[]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[0/3,0/3]-[0/0],esc,[1/3,0/3,0/3]-[]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[0/3,0/3]-[0/0],esc,[2/3,0/3,0/3]-[]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[0/3,0/3]-[0/0],esc,[3/3,0/3,0/3]-[]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[1/6,0/3]-[0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[2/6,0/3]-[0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[3/6,0/3]-[0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[4/6,0/3]-[0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[0/6]-[2/3,0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[1/6]-[2/3,0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[2/6]-[2/3,0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[3/6]-[2/3,0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[0/6]-[3/3,0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[1/6]-[3/3,0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[2/6]-[3/3,0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[3/6]-[3/3,0/0]], L = 6
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[1/6]-[0/0,0/0],esc,[0/3,1/6]-[0/0]], L = 6
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[1/6]-[0/0,0/0],esc,[1/3,1/6]-[0/0]], L = 6
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[1/6]-[0/0,0/0],esc,[2/3,1/6]-[0/0]], L = 6
%@ ;  Path = [sta,[1/3]-[0/0,0/0],sta,[1/6]-[0/0,0/0],esc,[3/3,1/6]-[0/0]], L = 6
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[1/6,0/3]-[0/0],esc,[0/3,...]-[]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[1/6,0/3]-[0/0],esc,[1/3,...]-[]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[1/6,0/3]-[0/0],esc,[2/3,...]-[]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[1/6,0/3]-[0/0],esc,[3/3,...]-[]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[2/6,0/3]-[0/0],des,[0/6]-[2/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[2/6,0/3]-[0/0],des,[1/6]-[2/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[2/6,0/3]-[0/0],des,[2/6]-[2/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[2/6,0/3]-[0/0],des,[3/6]-[2/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[3/6,0/3]-[0/0],des,[0/6]-[3/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[3/6,0/3]-[0/0],des,[1/6]-[3/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[3/6,0/3]-[0/0],des,[2/6]-[3/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[3/6,0/3]-[0/0],des,[3/6]-[3/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[4/6,0/3]-[0/0],des,[0/6]-[4/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[4/6,0/3]-[0/0],des,[1/6]-[4/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[4/6,0/3]-[0/0],des,[2/6]-[4/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[1/3,0/3]-[0/0],sta,[4/6,0/3]-[0/0],des,[3/6]-[4/6,...]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[0/6]-[2/3,0/0],esc,[2/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[0/6]-[2/3,0/0],esc,[3/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[0/6]-[2/3,0/0],esc,[4/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[0/6]-[2/3,0/0],esc,[5/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[1/6]-[2/3,0/0],esc,[2/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[1/6]-[2/3,0/0],esc,[3/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[1/6]-[2/3,0/0],esc,[4/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[2/3,0/3]-[0/0],des,[1/6]-[2/3,0/0],esc,[5/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[0/6]-[3/3,0/0],esc,[3/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[0/6]-[3/3,0/0],esc,[4/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[0/6]-[3/3,0/0],esc,[5/6,...]-[0/0]], L = 8
%@ ;  Path = [sta,[0/3]-[0/0,0/0],esc,[3/3,0/3]-[0/0],des,[0/6]-[3/3,0/0],esc,[6/6,...]-[0/0]], L = 8
%@ ;  ...
%@ caught: error(syntax_error(incomplete_reduction),read_term/3:1)

%?- length(Path, 2), phrase(path([0/0]-[0/0, 0/0]), Path), Path = [A,S].
%@    Path = [sta,[0/3]-[0/0,0/0]], A = sta, S = [0/3]-[0/0,0/0]
%@ ;  Path = [sta,[1/3]-[0/0,0/0]], A = sta, S = [1/3]-[0/0,0/0]
%@ ;  Path = [sta,[2/3]-[0/0,0/0]], A = sta, S = [2/3]-[0/0,0/0]
%@ ;  Path = [sta,[3/3]-[0/0,0/0]], A = sta, S = [3/3]-[0/0,0/0]
%@ ;  false.

%% Let's try a syntactically more elegant version of the above,
%% exploiting reified conjunction ','/3 from library(reif).

repath(debug(_)) --> []. % enable termination with debug info
repath(_) --> []. % a convenience for testing; path can stop at any time
repath(declare_mtd(_)) --> [].
repath(S0) --> { tpartition(state0_decision_regrettable(S0), [esc,sta,des], _, WontRegret),
		 (   WontRegret = [E|_], % NB: The ordering of [esc,sta,des] is *essential*
		     state0_decision_state(S0, E, S)
		 ;   WontRegret = [], % There aren't any non-regrettable decisions..
		     E = stop,        % ..other than STOP.
		     MTD = todo, % TODO: Actually obtain the MTD
		     S = declare_mtd(MTD)
		 )
	       }, [E, S],
	       repath(S).

%?- length(Path, 2), phrase(repath([0/0]-[0/0, 0/0]), Path), Path = [A,S].
%@    Path = [sta,[0/3]-[0/0,0/0]], A = sta, S = [0/3]-[0/0,0/0]
%@ ;  Path = [sta,[1/3]-[0/0,0/0]], A = sta, S = [1/3]-[0/0,0/0]
%@ ;  Path = [sta,[2/3]-[0/0,0/0]], A = sta, S = [2/3]-[0/0,0/0]
%@ ;  Path = [sta,[3/3]-[0/0,0/0]], A = sta, S = [3/3]-[0/0,0/0]
%@ ;  false.

%% TODO: Explore the goals inside cpath//2's {...}
%% The lesson here is that Prolog is not Python!
onestep(S0, E, S) :-
    if_(( % I believe this constructs a reified conjunction of reified conjunctions..
	       %(state0_decision_regrettable(S0, esc), E = esc), % G_esc
	       (state0_decision_regrettable(S0, sta), E = sta), % G_sta
	       (state0_decision_regrettable(S0, des), E = des)  % G_des
	   ),
	%% If all of [esc,sta,des] are regrettable, we STOP:
	(   E = stop,
	    MTD = todo, % TODO: Actually obtain the MTD
	    S = declare_mtd(MTD)
	),
	%% If any of the goals G_* above reified to FALSE,
	%% then I am hoping this branch will be entered
	%% with E instantiated to the first of [esc,sta,des]
	%% that is non-regrettable.
	(   must_be(atom, E),
	    state0_decision_state(S0, E, S)
	)
       ).

%?- onestep([0/0]-[0/0, 0/0], E, S).
%@ caught: error(instantiation_error,must_be/2)

%% I honestly think there's a good reason to introduce a 

%% TODO: Try to retain the clarity of repath//1 (contrast with path//1's
%%       spaghetti-nesting of multiple if_/3's!), but using reified conjunction,
%%       to exploit the efficiencies of early failures.
cpath(debug(_)) --> []. % enable termination with debug info
cpath(_) --> []. % a convenience for testing; path can stop at any time
cpath(declare_mtd(_)) --> [].
cpath(S0) --> { if_(( % I believe this constructs a reified conjunction of reified conjunctions..
			(E = esc, state0_decision_regrettable(S0, E)), % G_esc
			(E = sta, state0_decision_regrettable(S0, E)), % G_sta
			(E = des, state0_decision_regrettable(S0, E))  % G_des
		    ),
		    %% If all of [esc,sta,des] are regrettable, we STOP:
		    (	E = stop,
			MTD = todo, % TODO: Actually obtain the MTD
			S = declare_mtd(MTD)
		    ),
		    %% If any of the goals G_* above reified to FALSE,
		    %% then I am hoping this branch will be entered
		    %% with E instantiated to the first of [esc,sta,des]
		    %% that is non-regrettable.
		    state0_decision_state(S0, E, S)
		   ) % end of if/3
	       }, [E, S],
	       cpath(S).

%?- maplist(#=(3), [3,3]).
%@    true.
%@ false.

%?- length(Path, 2), phrase(cpath([0/0]-[0/0, 0/0]), Path), Path = [A,S].
%@    Path = [esc,[0/3,0/0]-[0/0]], A = esc, S = [0/3,0/0]-[0/0] % <- These
%@ ;  Path = [esc,[1/3,0/0]-[0/0]], A = esc, S = [1/3,0/0]-[0/0] % <- solutions
%@ ;  Path = [esc,[2/3,0/0]-[0/0]], A = esc, S = [2/3,0/0]-[0/0] % <- are
%@ ;  Path = [esc,[3/3,0/0]-[0/0]], A = esc, S = [3/3,0/0]-[0/0] % <- wrong.
%@ ;  Path = [sta,[0/3]-[0/0,0/0]], A = sta, S = [0/3]-[0/0,0/0]
%@ ;  Path = [sta,[1/3]-[0/0,0/0]], A = sta, S = [1/3]-[0/0,0/0]
%@ ;  Path = [sta,[2/3]-[0/0,0/0]], A = sta, S = [2/3]-[0/0,0/0]
%@ ;  Path = [sta,[3/3]-[0/0,0/0]], A = sta, S = [3/3]-[0/0,0/0]
%@ ;  false.

%% How did the above 'pass', given that esc from [0/0]-[0/0,0/0] is truly regrettable?
%?- state0_decision_regrettable([0/0]-[0/0,0/0], esc, T).
%@    T = true
%@ ;  false.

%% But shouldn't I have anticipated possible regret upon escalating from 0/0?
%?- regret(esc, [1/3,0/0]).
%@    true
%@ ;  false.

%?- N in 1..8, indomain(N), length(Path, N), phrase(path([]-[0/0, 0/0]), Path).
%@    N = 2, Path = [esc,[0/3]-[0/0]]
%@ ;  N = 2, Path = [esc,[1/3]-[0/0]]
%@ ;  N = 2, Path = [esc,[2/3]-[0/0]]
%@ ;  N = 2, Path = [esc,[3/3]-[0/0]]
%@ ;  N = 4, Path = [esc,[0/3]-[0/0],esc,[0/3,0/3]-[]]
%@ ;  N = 4, Path = [esc,[0/3]-[0/0],esc,[1/3,0/3]-[]]
%@ ;  N = 4, Path = [esc,[0/3]-[0/0],esc,[2/3,0/3]-[]]
%@ ;  N = 4, Path = [esc,[0/3]-[0/0],esc,[3/3,0/3]-[]]
%@ ;  N = 4, Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0]]
%@ ;  N = 4, Path = [esc,[1/3]-[0/0],sta,[2/6]-[0/0]]
%@ ;  N = 4, Path = [esc,[1/3]-[0/0],sta,[3/6]-[0/0]]
%@ ;  N = 4, Path = [esc,[1/3]-[0/0],sta,[4/6]-[0/0]]
%@ ;  N = 6, Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[]]
%@ ;  N = 6, Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[]]
%@ ;  N = 6, Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[2/3,1/6]-[]]
%@ ;  N = 6, Path = [esc,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[3/3,1/6]-[]]
%@ ;  false.

%% I wonder if the search, failing to find a path with length L,
%% searches for one with length L+1, then L+2, ad infinitum.
%% If this is the case, does that vitiate my goal, rather than
%% the underlying mechanics? I think so!
%% Then I need a way to query these designs without such
%% explosions.

%?- phrase(path([3/3,1/6]-[]), Path).
%@    Path = []
%@ ;  false.
