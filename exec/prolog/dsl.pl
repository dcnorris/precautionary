/*

On the way toward a domain-specific language (DSL) for dose-escalation designs,
we implement the rules of the classical 3+3 design [1,2], formulated in term of
*prohibitions* on certain future events considered to be 'regrettable'.

1. Korn EL, Midthune D, Chen TT, Rubinstein LV, Christian MC, Simon RM.
   A comparison of two phase I trial designs. Stat Med. 1994;13(18):1799-1806.
   doi:10.1002/sim.4780131802

2. Norris DC. What Were They Thinking? Pharmacologic priors implicit in a choice
   of 3+3 dose-escalation design. arXiv:201205301 [statME]. December 24, 2020.
   https://arxiv.org/abs/2012.05301

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

With a suitable interpreter defined, it might be possible generally to operate
in a mode of https://en.wikipedia.org/wiki/Inductive_programming, to abstract
approximating 'regret' rules from the paths of any dose-escalation design.

If this proves possible, then the whole class of dose-escalation trial designs
would be effectively subsumed under this regret-constraint framework.
This, in turn, will be useful if this framework confers a *generalizability*
upon these designs (and upon their analysis) similar to that enjoyed by other
constraint-type formulations throughout mathematics and the sciences. Thus,
one may hope that dose-escalation designs defined in regret-constraint terms
will readily generalize to rolling enrollment [3,4] -- and perhaps even to
dose-TITRATION as in the 3+3/PC design [5].

3. Skolnik JM, Barrett JS, Jayaraman B, Patel D, Adamson PC. Shortening the
   timeline of pediatric phase I trials: the rolling six design. J Clin Oncol.
   2008;26(2):190-195. doi:10.1200/JCO.2007.12.7712

4. Frankel PH, Chung V, Tuscano J, et al. Model of a Queuing Approach for
   Patient Accrual in Phase 1 Oncology Studies. JAMA Network Open.
   2020;3(5):e204787-e204787. doi:10.1001/jamanetworkopen.2020.4787

5. Norris DC. Precautionary Coherence Unravels Dose Escalation Designs.
   bioRxiv. December 29, 2017. doi:10.1101/240846

*/

%% What follows is a demonstration that a regret-constraint formulation can
%% recapitulate  the common 3+3 design. Apart from using standard libraries,
%% this program is self-contained.
:- use_module(library(lists)).
:- use_module(library(clpz)).
:- use_module(library(reif)).
:- use_module(library(si)).
:- use_module(library(dcgs)).
:- use_module(library(error)).
:- use_module(library(lambda)).

%% During the course of a dose-escalation trial, at any given dose-level
%% there is some current tally T/N, T,N ∈ ℕ recording T toxic responses
%% ('toxicities') out of N participants enrolled at that dose. Enrollment
%% of a further 'cohort' at this dose-level increases the denominator and
%% possibly also the numerator. Initially, we constrain enrollment to
%% cohorts of size 3, and limit total enrollment at any one dose to 6.
%% (The 'toxicity assessment period' that must elapse before a toxicity
%% determination can be made is NOT explicitly modeled, and is thus elided
%% as if it were instantaneous.)
enroll(T0/N0, T1/N1) :-
    #Nnew #= 3, % hard-coding cohorts of 3 for now
    #N1 #= #N0 + #Nnew,
    Tnew in 0..Nnew, indomain(Tnew),
    #T1 #= #T0 + #Tnew,
    #N1 #=< 6. % Maximum enrollment per cohort
    
%?- enroll(0/0, T/N).
%@    T = 0, N = 3
%@ ;  T = 1, N = 3
%@ ;  T = 2, N = 3
%@ ;  T = 3, N = 3.

%% We model the state of the whole trial (comprising multiple dose-levels)
%% as a pair of lists of tallies, the left-hand list in *descending* dose
%% order with the 'current dose' at its head, and the right-hand list in
%% *ascending* order, with the next-higher dose at its head. Thus doses
%% adjacent to the 'current dose' are immediately accessible. Dose-skipping
%% is generally not done in dose-escalation trials, in which case there are
%% 3 main dose-escalation decisions: ESCalate, STAy and DeEScalate:
state0_decision_state(Ls - [R0|Rs], esc, [R|Ls] - Rs) :- enroll(R0, R).
state0_decision_state([L0|Ls] - Rs, sta, [L|Ls] - Rs) :- enroll(L0, L).
state0_decision_state([L,D0|Ls] - Rs, des, [D|Ls] - [L|Rs]) :- enroll(D0, D).

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

%% Note how a trial starts naturally from the obvious 'initial state':
%?- state0_decision_state([] - [0/0, 0/0], esc, S).
%@    S = [0/3]-[0/0]
%@ ;  S = [1/3]-[0/0]
%@ ;  S = [2/3]-[0/0]
%@ ;  S = [3/3]-[0/0]
%@ ;  false.

%% ~~~~~~~~~~ Abstracting the 'regrets' implicit in 3+3 ~~~~~~~~~~

%% regret(A, [Q|Qs]) holds if we regret having taken action A if it
%% results in a tally history [Q|Qs] -- which reads backwards in time.
%% This allows for the relevant history to go back either 1 or 2 tallies,
%% which is fully suffient for practical purposes. (Regrets based on long
%% histories would hardly be intuitive for clinical trialists, and a key
%% aim here is of course to capture the definitive intuitions.)

%% This says that, regardless of the resulting tally, we would regret
%% having de-escalated from an 'insufficiently toxic' tally of 1-/3+.
regret(des, [_,T/N]) :-
    #T #=< 1 #/\ #N #>= 3.
%% TODO: Try to tighten this up, so that we regret 0/3 toxicities after
%%       having de-escalated from a 1/3.
%% TODO: Might this regret already be expressed via the preference for
%%       'sta' when possible? Or do we need to regret 'des' in order to
%%       trigger the 'stop' decision?
%% NOTE: Some of these decisions may best be deferred until our attempts
%%       at generalization can offer selective principles.

%% We regret having stayed at the current dose upon seeing a 5+/6 tally.
regret(sta, [T/N|_]) :-
    #N #= 6 #/\ #T #>= 5. % regret excess toxicity

%% It is in regretting ESCALATION decisions that the clinical intuition
%% expresses itself most strongly. The core intuition here is almost
%% 'medico-legal' in character. What is regretted is the occurrence of
%% (certain degreees of) toxicity without 'sufficient excuse' provided
%% by experiences in the next-lower dose level. (Cf. 'safe harbor'.)
%% Thus, for example, we regret ANY amount of toxicity at dose D+1
%% if we lack a basis of 0/3 or no more than 1/6 toxicities at dose D.

regret(esc, [T/N, T0/N0]) :-
    #N #>= #T, % condition for T/N to be a valid tally
    (	#T #> 0, % We will regret even 1 toxicity when escalating after..
	 #N0 #< 3 % having enrolled less than 3 at previous dose.
    ;	#T #> 1, % We regret >1 toxicities after..
	#T0 * 6 #> N0 % having seen tox rate T0/N0 > 1/6 at prev dose.
    ;	#T #>= 5 % We regret net 5+ toxicities after any escalation.
	%% NB:
	%% The need for this final case is not obvious. It arises when
	%% we have accumulated 0/6 toxicities at current dose, having
	%% de-escalated from a higher dose with tally 2+/3. That this
	%% case is so nontrivial perhaps points toward a need for some
	%% new concept to render it 'obvious' or 'natural'.
    ).

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

%% I believe the semantics of 'regret' are untenable at the level of if_/3
%% application. What I need is a goal that reifies to 'true' for precisely
%% those decisions which are BOTH feasible AND non-regrettable, and otherwise
%% reifies to 'false'.
state0_decision_noregrets(S0, A, Truth) :-
    state_si(S0),
    regrettable_decision(A),
    if_(state0_decision_feasible(S0, A),
	(   state0_decision_state(S0, A, S),
	    S0 = [T0/N0|_] - _, % TODO: Factor this pattern matching
	    S  = [T /N |_] - _, %       into a regret/3 predicate?
	    %% I introduce the (->) below to avert backtracking over
	    %% possibly multiple scenarios for regret -- one is enough!
    	    regret(A, [T/N, T0/N0]) -> Truth = false
	;   Truth = true
	),
	Truth = false
       ).

%% These are the 3 dose-escalation decisions
%% to which the 'regret' concept applies:
regrettable_decision(esc).
regrettable_decision(sta).
regrettable_decision(des).

state0_decision_feasible(S0, A, Truth) :-
    state_si(S0),
    (	state0_decision_state(S0, A, _) -> Truth = true
    ;	Truth = false
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

%?- state0_decision_noregrets([1/3]-[0/0], E, T).
%@    E = esc, T = false
%@ ;  E = sta, T = true
%@ ;  false.

%% This code reflects the primacy of REGRET as the crucial user-level
%% concept shaping these designs.

cascading_decision_otherwise([], Os, _, _, _) :- maplist(call, Os).
cascading_decision_otherwise([D|Ds], Os, E, S0, S) :-
        if_(state0_decision_noregrets(S0, D),
            (   E = D,
                state0_decision_state(S0, E, S)
            ),
            cascading_decision_otherwise(Ds, Os, E, S0, S)).

%path(_) --> []. % a convenience for testing; path can stop at any time
path(declare_mtd(_)) --> [].
path(S0) --> { cascading_decision_otherwise([esc,sta,des],
                                            [E = stop, % ..and otherwise, STOP.
                                             MTD = todo, % TODO: Actually obtain MTD as integer >= 0.
                                             S = declare_mtd(MTD)], E, S0, S) },
	     [E, S],
	     path(S). % TODO: Implement declare_mtd possibility

%?- Path = [sta,[3/3]-[0/0],stop,declare_mtd(todo)], phrase(path([0/0]-[0/0]), Path).
%@    Path = [sta,[3/3]-[0/0],stop,declare_mtd(todo)]
%@ ;  false.

%% Right away, we can see that 'des' occurs when we should declare_mtd/1.
%% This makes me wonder whether letting 'des' be some kind of default
%% is wrong. Perhaps I should also have a concept of regretting 'des'?
%% Or else maybe there should be some *more active* effort to declare_mtd/1
%% ASAP.
%% I really do like the idea of defining the MTD implicitly through termination
%% of the dose-escalation, by 'running out of options'.
%% Does retaining this scheme require defining regret even for 'des'?
%% Or do I need to take special account of enrollment limits?
%% REMEMBER: I am hoping to make the definition of these trials depend on
%% only a very small set of rules. If I have to start dealing separately
%% with this-or-that type of regret, does that undermine this key aim?

%?- state0_decision_noregrets([0/3,0/3]-[], E, Truth).
%@    E = esc, Truth = false
%@ ;  E = sta, Truth = true
%@ ;  E = des, Truth = true
%@ ;  false.

%?- state0_decision_noregrets([0/6,0/3]-[], des, Truth).
%@    Truth = false
%@ ;  false.

%?- phrase(path([0/0]-[0/0]), Path).
%@    Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[0/6,0/3]-[],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[1/6,0/3]-[],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[0/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[1/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[2/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[3/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[0/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[1/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[2/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[3/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[1/6,0/3]-[],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[0/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[1/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[2/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[3/6]-[2/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[0/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[1/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[2/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[3/6]-[3/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[0/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[1/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[2/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[3/6]-[4/6],stop,declare_mtd(...)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[0/6]-[2/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[1/6]-[2/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[2/6]-[2/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[3/6]-[2/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[0/6]-[3/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[1/6]-[3/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[2/6]-[3/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[3/6]-[3/3],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[0/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[1/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[2/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[1/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[2/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[3/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[4/6,...]-[],stop,declare_mtd(...)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[2/3,1/6]-[],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[3/3,1/6]-[],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[2/6]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[3/6]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[4/6]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[2/3]-[0/0],stop,declare_mtd(todo)]
%@ ;  Path = [sta,[3/3]-[0/0],stop,declare_mtd(todo)]
%@ ;  false. % J=46 paths!

%?- J+\(length(D,1), maplist(=(0/0), D), findall(Path, phrase(path([]-D), Path), Paths), length(Paths, J)).
%@    J = 10.
%?- J+\(length(D,2), maplist(=(0/0), D), findall(Path, phrase(path([]-D), Path), Paths), length(Paths, J)).
%@    J = 46.
%?- J+\(length(D,3), maplist(=(0/0), D), findall(Path, phrase(path([]-D), Path), Paths), length(Paths, J)).
%@    J = 154.
%?- J+\(length(D,4), maplist(=(0/0), D), findall(Path, phrase(path([]-D), Path), Paths), length(Paths, J)).
%@    J = 442.
%?- J+\(length(D,5), maplist(=(0/0), D), findall(Path, phrase(path([]-D), Path), Paths), length(Paths, J)).
%@    J = 1162.
%?- J+\(length(D,6), maplist(=(0/0), D), findall(Path, phrase(path([]-D), Path), Paths), length(Paths, J)).
%@    J = 2890.
%% Yep!

/*

PLAN of further work ...

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

