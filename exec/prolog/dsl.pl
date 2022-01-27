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
:- use_module(library(time)).
:- use_module(library(debug)).

:- initialization(assertz(clpz:monotonic)).

%% During the course of a dose-escalation trial, at any given dose-level
%% there is some current tally T/N, T,N ∈ ℕ recording T toxic responses
%% ('toxicities') out of N participants enrolled at that dose. Enrollment
%% of a further 'cohort' at this dose-level increases the denominator and
%% possibly also the numerator. This predicate allows for C to be specified
%% as a (possibly singleton) range. Initially, we hard-code a limit of 6
%% on total enrollment at any one dose.
%% (The 'toxicity assessment period' that must elapse before a toxicity
%% determination can be made is NOT explicitly modeled, and is thus elided
%% as if it were instantaneous.)
enroll(C/M, T0/N0, T1/N1) :-
    Nnew in C, indomain(Nnew), % C is a (possibly singleton) integer range
    #N1 #= #N0 + #Nnew,
    Tnew in 0..Nnew, indomain(Tnew),
    #T1 #= #T0 + #Tnew,
    #N1 #=< M. % Maximum enrollment per cohort
    
%?- enroll(3/6, 0/0, T/N).
%@    T = 0, N = 3
%@ ;  T = 1, N = 3
%@ ;  T = 2, N = 3
%@ ;  T = 3, N = 3
%@ ;  false. % Why extra choice-point when C is a singleton?

%?- enroll((2..3)/6, 0/0, T/N).
%@    T = 0, N = 2
%@ ;  T = 1, N = 2
%@ ;  T = 2, N = 2
%@ ;  T = 0, N = 3
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
state0_decision_state(Ls - [H0|Hs], esc, [H|Ls] - Hs) :-
    cohort(C), maxenr(M), enroll(C/M, H0, H).
state0_decision_state([L0|Ls] - Hs, sta, [L|Ls] - Hs) :-
    cohort(C), maxenr(M), enroll(C/M, L0, L).
state0_decision_state([L,D0|Ls] - Hs, des, [D|Ls] - [L|Hs]) :-
    cohort(C), maxenr(M), enroll(C/M, D0, D).
:- discontiguous(cohort/1). % (to permit later additions)
cohort(3). % For the present application we use cohorts of size 3,
maxenr(6). % and enroll at most 6 patients at any one dose level.

%?- state0_decision_state([0/3, 1/6] - [0/0, 0/0], esc, Ls - Hs).
%@    Ls = [0/3,0/3,1/6], Hs = [0/0]
%@ ;  Ls = [1/3,0/3,1/6], Hs = [0/0]
%@ ;  Ls = [2/3,0/3,1/6], Hs = [0/0]
%@ ;  Ls = [3/3,0/3,1/6], Hs = [0/0]
%@ ;  false.

%?- state0_decision_state([1/3, 1/6] - [0/0, 0/0], sta, Ls - Hs).
%@    Ls = [1/6,1/6], Hs = [0/0,0/0]
%@ ;  Ls = [2/6,1/6], Hs = [0/0,0/0]
%@ ;  Ls = [3/6,1/6], Hs = [0/0,0/0]
%@ ;  Ls = [4/6,1/6], Hs = [0/0,0/0]
%@ ;  false.

%?- state0_decision_state([2/3, 0/3] - [0/0, 0/0], des, Ls - Hs).
%@    Ls = [0/6], Hs = [2/3,0/0,0/0]
%@ ;  Ls = [1/6], Hs = [2/3,0/0,0/0]
%@ ;  Ls = [2/6], Hs = [2/3,0/0,0/0]
%@ ;  Ls = [3/6], Hs = [2/3,0/0,0/0]
%@ ;  false.

%% Note how a trial starts naturally from the obvious 'initial state':
%?- S+\(Zero = []-[0/0,0/0], state0_decision_state(Zero, esc, S)).
%@    S = [0/3]-[0/0]
%@ ;  S = [1/3]-[0/0]
%@ ;  S = [2/3]-[0/0]
%@ ;  S = [3/3]-[0/0]
%@ ;  false.

%% Note the manner in which Ls-Hs maps to a dose-ordered tally list:
state_tallies(Ls-Hs, Qs) :- reverse(Ls,Js), append(Js,Hs,Qs).

%?- state_tallies([]-[1/6,2/3], Qs).
%@    Qs = [1/6,2/3].

%% ~~~~~~~~~~ Abstracting the 'regrets' implicit in 3+3 ~~~~~~~~~~

%% regret(A, [Q|Qs]) holds if we regret having taken action A if it
%% results in a tally history [Q|Qs] -- which reads backwards in time.
%% This allows for the relevant history to go back either 1 or 2 tallies,
%% which is fully suffient for practical purposes. (Regrets based on long
%% histories would hardly be intuitive for clinical trialists, and a key
%% aim here is of course to capture the definitive intuitions.)

:- discontiguous(regret/2).
%% This says that, regardless of the resulting tally, we would regret
%% having de-escalated from an 'insufficiently toxic' tally of 1-/3+.
regret(des, [T/N, T0/N0]) :-
    #T0 #=< 1 #/\ #N0 #>= 3, % prev dose no more than moderately toxic
    #N #> 0, #T * 6 #> #N.   % new toxicity rate is lower than 1/6
%% TODO: Try to tighten this up, so that we regret 0/3 toxicities after
%%       having de-escalated from a 1/3.
%% TODO: Might this regret already be expressed via the preference for
%%       'sta' when possible? Or do we need to regret 'des' in order to
%%       trigger the 'stop' decision?
%% NOTE: Some of these decisions may best be deferred until our attempts
%%       at generalization can offer selective principles.

%% It is in regretting ESCALATION decisions that the clinical intuition
%% expresses itself most strongly. The core intuition here is almost
%% 'medico-legal' in character. What is regretted is the occurrence of
%% (certain degreees of) toxicity without 'sufficient excuse' provided
%% by experiences in the next-lower dose level. (Cf. 'safe harbor'.)
%% Thus, for example, we regret ANY amount of toxicity at dose D+1
%% if we lack a basis of 0/3 or no more than 1/6 toxicities at dose D.

regret(esc, [T/N, T0/N0]) :- % regret escalating without `safe harbor'
    #N #>= #T, % condition for T/N to be a valid tally
    (   #T #> 0, % We will regret even 1 toxicity when escalating after..
        #N0 #< 3 % having enrolled less than 3 at previous dose.
    ;   #T #> 1, % We regret >1 toxicities after..
        #T0 * 6 #> N0 % having seen tox rate T0/N0 > 1/6 at prev dose.
    ).

%% Finally, we regret ANY decision that results in a 5+/_ tally:
regret(_, [T/N|_]) :-
    #T #>= 5 #/\ #N #>= #T. % regret excess toxicity

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
state0_decision_noregrets(S0, E, Truth) :-
    state_si(S0),
    regrettable_decision(E),
    if_(state0_decision_feasible(S0, E),
        (   state0_decision_state(S0, E, S),
            S0 = [T0/N0|_] - _, % TODO: Factor this pattern matching
            S  = [T /N |_] - _, %       into a regret/3 predicate?
            %% The (->) below averts backtracking over possibly multiple
            %% scenarios for regret. Regret being idempotent, 1 suffices.
            regret(E, [T/N, T0/N0]) -> Truth = false
        ;   Truth = true
        ),
        Truth = false
       ).

regrettable_decision(esc). % These are the 3 dose-escalation decisions
regrettable_decision(sta). % to which the `regret' concept applies.
regrettable_decision(des). % Note in particular the absence of `stop'.

%% TODO: Does the (->)/2 below really engender distrust,
%%       even though we can invoke with E uninstantiated?
state0_decision_feasible(S0, E, Truth) :-
    state_si(S0),
    regrettable_decision(E), % guarantees E is ground
    (   state0_decision_state(S0, E, _) -> Truth = true
    ;   Truth = false
    ).

%% TODO: State the monotonic execution constraint!
tally_si(T/N) :-
    maxenr(MaxN),
    N in 0..MaxN, indomain(N), % TODO: No indomain/1 needed in an _si predicate.
    T in 0..N,
    indomain(T). % as an adapted process, state0_decision_noregrets/3 demands ground state

%% To explore fully how the natural-language protocol descriptions relate
%% to our formal trial definition, we will require a state_si/1 predicate
%% capable of generating fairly enumerated states.
state_si(Ls - Hs) :-
    D in 1..7,     % TODO: A catastrophic 'wart' on otherwise elegant code?
    indomain(D),
    Dcur in 0..D,  %       Does all this merely remind us how indispensable
    indomain(Dcur),
    #Dhi #= D - Dcur,
    length(Ls, Dcur), %     clpz is for logic programming? Or does it
    length(Hs, Dhi),  %     suggest we need a smarter append/3?
    %% ~~ Everything below this point seemed so delightfully elegant ~~
    %%%%length(Doses, D), % we prespecify a number D of dose levels
    append(Ls, Hs, Doses), % the 'current dose' splits Doses somewhere
    maplist(tally_si, Doses). % all doses bear a tally

%?- state0_decision_noregrets([1/3]-[0/0], E, T).
%@    E = esc, T = false
%@ ;  E = sta, T = true
%@ ;  E = des, T = false
%@ ;  false.

%% Now we can ask rather general queries based on partially instantiated states!
%% The following says what we can do without regret:
%?- setof(E, Q^state0_decision_noregrets([1/3,Q]-[0/0], E, true), Okay).
%@    Okay = [sta].

%?- setof(E, Q^(regrettable_decision(E), state0_decision_noregrets([1/3,Q]-[0/0], E, true)), Okay).
%@    Okay = [sta].

%?- setof(E, T^(T in 0..6, indomain(T), state0_decision_noregrets([T/6,0/0]-[0/0], E, true)), Okay).
%@    Okay = [des,esc].

%?- setof(E, T^(T in 0..6, state0_decision_noregrets([T/6,0/0]-[0/0], E, true)), Okay).
%@ false.

%% Other demonstrations we can effect include:
%?- time(setof(E, P^R^state0_decision_noregrets([0/3,P]-[R], E, true), Okay)).
%@    % CPU time: 157.494s
%@    Okay = [esc,sta].

%?- time(setof(T:E, P^R^state0_decision_noregrets([T/6,P]-[R], E, true), Okay)).
%@    % CPU time: 1008.281s
%@    Okay = [0:esc,1:esc,2:des,3:des,4:des,5:des,6:des].

%% 1. The first action is 'escalating' from a zero dose!
%?- length(Levels, 3), maplist(=(0/0), Levels), S0 = []-Levels, setof(E, state0_decision_noregrets(S0, E, true), Begin).
%@    Levels = [0/0,0/0,0/0], S0 = []-[0/0,0/0,0/0], Begin = [esc].

%% 2. If 0/3 DLTs, escalate
%?- setof(E, P^R^state0_decision_noregrets([0/3,P]-[R], E, true), Okay).
%@    Okay = [esc,sta]. % NB: preference for esc > sta would bind here

%% 3. If 1/3 DLTs, treat an additional 3 at current dose level
%?- setof(E, P^R^state0_decision_noregrets([1/3,P]-[R], E, true), Okay).
%@    Okay = [sta].

%% 4. If 1/6 DLTs, escalate
%?- setof(E, P^R^state0_decision_noregrets([1/6,P]-[R], E, true), Okay).
%@    Okay = [esc].

%% 4b. Edge case omitted: we are at top dose
%?- state0_decision_noregrets([1/6,_]-[], E, true).
%@ false. % Cannot find a dose-escalation decision we won't regret ==> STOP

%% TODO: May have to deal separately with cases where previous dose
%%       had either 3 or 6 patients treated.
%% 5. If T >= 2, then MTD has been exceeded
%?- time(setof(E, N^P^R^(#T #>= 2, state0_decision_noregrets([T/N,P]-[R], E, true)), Okay)).
%@    % CPU time: 2145.148s
%@    Okay = [des], T = 2
%@ ;  % CPU time: 0.001s
%@    Okay = [des], T = 3
%@ ;  % CPU time: 0.000s
%@    Okay = [des], T = 4
%@ ;  % CPU time: 0.000s
%@    Okay = [des], T = 5
%@ ;  % CPU time: 0.000s
%@    Okay = [des], T = 6.

recommendation_exceeds_mtd(NDoses) :-
    #NDoses #> 0, % Number of doses must be a positive integer. 
    length(Ds,NDoses), maplist(=(0/0), Ds), % From initial state D,
    phrase(path([]-Ds), Path),  % .. we seek a Path of the trial
    phrase((..., [Ls-_], ...,   % .. on which state Ls-_ appeared,
	    [recommend_dose(Rec)] %  and the recommended dose Rec
	   ), Path),
    length(Ls,X), Rec #>= X, % .. was at least the current dose X,
    Ls = [T/_|_], #T #> 1.   % .. yet X `exceeded MTD' per protocol.

%?- time((N in 1..3, indomain(N), recommendation_exceeds_mtd(N))).
%@    % CPU time: 20.626s
%@ false.
%@    % CPU time: 0.003s
%@ caught: error(existence_error(procedure,report_excessive_mtd/1),report_excessive_mtd/1)
%@    % CPU time: 0.001s
%@ caught: error(existence_error(procedure,report_excessive_mtd/1),report_excessive_mtd/1)
%@    % CPU time: 1195.630s % N in 1..7
%@ false.
%@    % CPU time: 480.133s % N in 1..6
%@ false.

%% To support the paper, a free query (not RHS of Horn clause) serves best:
badness :-
    D in 1..7, indomain(D), % For trials of practical size (up to 7 doses)
    InitN = []-Ds, length(Ds, D), % .. that start from the lowest dose
    maplist(=(0/0), Ds),          % .. with no prior toxicity information,
    phrase(path(InitN), Path),    % can we find a propery-violating Path
    phrase((..., [Ls-_], ...,     % .. on which a state Ls-_ appears,
            [recommend_dose(Rec)] % .. such that the recommended dose Rec
           ), Path),
    length(Ls,X), Rec #>= X,      % .. was at least the current dose X,
    Ls = [T/_|_], #T #> 1.        % .. yet X `exceeded MTD' per protocol?

%?- time(badness).
%@    % CPU time: 1262.907s
%@ false.

%% TODO: Verify the 'definition' of MTD.
rec_not_as_defined :-
    D in 1..7, indomain(D), % For trials of usual size (up to D=7 doses)
    InitD = []-Ds, length(Ds, D), % .. that start from the lowest dose
    maplist(=(0/0), Ds),          % .. with no prior toxicity information,
    phrase(path(InitD), Path),    % does any Path exist
    phrase((..., [S, stop,        % .. on which a terminal state S appears,
                  recommend_dose(Rec)] % .. and recommended dose is Rec
                 ), Path),
    state_tallies(S, Qs), % such that the dose-ordered tally list ..
    length(OKs, Rec), append(OKs, TOXs, Qs), % split on the recommendation
    (   nth0(Rec, [nil|Qs], T/N), % has the recommended dose
        #T #> 1 #\/ #N #< 6 % .. being too toxic or insufficiently tried
    ;   \+maplist(\Q^(Q=T/N, #T #> 1 #\/ #N #< 6), % .. OR the same fails
                  TOXs) % .. to hold for any of the above-Rec doses?
    ).

%?- time(rec_not_as_defined).
%@    % CPU time: 1359.726s
%@ false.

%% ALL DIRECTIONS
%% Looking ONE STEP FORWARD ..
%% NB: This required reverting to 'master'; see #1246
%?- setof(E, Path^(phrase(path([]-[0/0,0/0,0/0]), Path), phrase((...,[[0/3,0/3,0/3]-[], E], ...), Path)), Es).
%@    Es = [sta].

%% Looking BACKWARD ..
%?- setof(E0, Path^(phrase(path([]-[0/0,0/0,0/0]), Path), phrase((...,[E0, [0/3,0/3,0/3]-[]],...), Path)), E0s).
%@    E0s = [esc].

%% Looking FAR FORWARD ..
%?- setof(Rec, Path^(phrase(path([]-[0/0,0/0,0/0]), Path), phrase((...,[[2/6,0/3,0/3]-[]],...,[recommend_dose(Rec)]), Path)), Recs).
%@    Recs = [0,1,2].

%% A fresh way to see PROTOCOL VIOLATION:
%?- phrase(path([]-[0/0,0/0,0/0,0/0]), Path), phrase((..., [[0/3,0/3,0/3]-[0/0], _, [2/6,0/3,0/3]-[0/0]], ...), Path).
%@ false.

%% Finally, COMPLETE ENUMERATION of all trial paths:
%?- J+\time((D = 7, length(Ds, D), maplist(=(0/0), Ds), setof(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J))).
%@    % CPU time: 731.888s % D=7 (rebis-dev, i7-4790 @ 4.00GHz)
%@    J = 6922.
%@    % CPU time: 304.874s % D=6 (rebis-dev)
%@    J = 2890.
%@    % CPU time: 120.060s % D=5 (rebis-dev)
%@    J = 1162.
%@    % CPU time: 43.237s  % D=4 (rebis-dev)
%@    J = 442.
%@    % CPU time: 222.425s % D=4 (master)
%@    J = 442.

%% LIVENESS & OTHER PROPERTIES

%% The trial always finishes with a recommendation between 0 and D.
recommends_wrongly(D) :-
    #D #> 0,
    InitD = []-Ds, length(Ds, D), maplist(=(0/0), Ds),
    phrase(path(InitD), Path), % Does any Path exist
    (   phrase((...,           % .. which ever
                [recommend_dose(Rec)], % recommends a dose
                [_], ...), Path) % .. and then continues?
    ;   \+ phrase((..., % .. or fails to end with a rec?
                   [recommend_dose(Rec)]), Path)
    ;   phrase((..., % .. or ever recommends a dose
                [recommended_dose(Rec)], ...), Path),
        #\ (Rec in 0..D) % .. that is not an integer from 0 to D?
    ).

%?- time((D in 1..4, indomain(D), recommends_wrongly(D))).
%@    % CPU time: 62.902s
%@ false.

nondeterminism(D) :-
    #D #> 0,
    length(Ds, D), maplist(=(0/0), Ds),
    phrase(path([]-Ds), Path),
    phrase((..., [Ls-Hs, E1], ...), Path),
    phrase((..., [Ls-Hs, E2], ...), Path),
    E1 \= E2.

%?- time((D in 1..3, indomain(D), nondeterminism(D))).
%@    % CPU time: 19.792s
%@ false.
%@ false.

%?- cohort(C).
%@    C = 1
%@ ;  C = 2
%@ ;  C = 3.

%?- state0_decision_state([0/3,0/3,0/3]-[], sta, [T/N,0/3,0/3]-[]).
%@    T = 0, N = 4
%@ ;  T = 1, N = 4
%@ ;  T = 0, N = 5
%@ ;  T = 1, N = 5
%@ ;  T = 2, N = 5
%@ ;  T = 0, N = 6
%@ ;  T = 1, N = 6
%@ ;  T = 2, N = 6
%@ ;  T = 3, N = 6
%@ ;  false.

%% ROLLING ENROLLMENT
regret(sta, [_, T/_]) :- % NB: a backward-looking regret
    #T #>= 2. % regret staying at a too-toxic dose we'd never recommend
regret(esc, [3/N|_]) :- N in 4..6, indomain(N).

%% TODO: Show that neither of these additional regrets
%%       eliminates any path of the trial:
%?- J+\(D = 7, length(Ds, D), maplist(=(0/0), Ds), setof(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J)).
%@    J = 6922. % no paths have been lost

:- discontiguous(cohort/1).
cohort(1).
cohort(2).
cohort(3).

%% More flexible enrollment as above now enables de-escalation
%% before dosing participant C in Figure 1 of the paper:
%?- phrase(path([0/3,0/3,0/3]-[]), [sta, [2/5,0/3,0/3]-[], des, [0/6,0/3]-[2/5], stop, recommend_dose(2)]).
%@    true
%@ ;  ...

%?- phrase(path([0/3,0/3,0/3]-[]), Path).
%@    Path = [sta,[0/5,0/3,0/3]-[],stop,recommend_dose(3)]
%@ ;  Path = [sta,[1/5,0/3,0/3]-[],stop,recommend_dose(2)]
%@ ;  Path = [sta,[2/5,0/3,0/3]-[],des,[0/5,0/3]-[2/5],stop,recommend_dose(2)]
%@ ;  Path = [sta,[2/5,0/3,0/3]-[],des,[1/5,0/3]-[2/5],stop,recommend_dose(1)]
%@ ;  Path = [sta,[2/5,0/3,0/3]-[],des,[2/5,0/3]-[2/5],des,[0/5]-[2/5,2/5],stop,recommend_dose(1)]
%@ ;  ... 

%% CHALLENGE: Can I write a query that regurgitates the entirety of Korn '94 protocol?
%?- setof(Q:E, P^R^state0_decision_noregrets([Q,P]-[0/0], E, true), Okay).
%@    Okay = [0/0:des,0/0:sta,0/1:des,0/1:sta,0/2:des,0/2:sta,0/3:esc,0/3:sta,0/4:esc,... : ...|...].

%?- time(setof(Q, N^P^R^(Q = T/N, N in 0..6, T in 0..6, T #=< N, state0_decision_noregrets([Q,P]-[R], esc, true)), Okay)).
%@    % CPU time: 1320.068s
%@    Okay = [0/3,0/4,0/5,0/6], T = 0
%@ ;  % CPU time: 0.000s
%@    Okay = [1/6], T = 1.

%?- T in 0..6, indomain(T), state0_decision_noregrets([T/N,0/0]-[0/0], E, true).

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
path(recommend_dose(_)) --> [].
path(S0) --> { cascading_decision_otherwise([esc,sta,des],
                                            [E = stop,
                                             stopstate_rec(S0, Rec),
                                             S = recommend_dose(Rec)], E, S0, S) },
	     [E, S],
	     path(S).

%% In order to validate the 'global' constraints implicit in Korn'94,
%% we must at last define 'the' MTD!
%% TODO: Refine and condense this, as much as possible.
stopstate_rec(S, Rec) :-
    state_si(S),
    (   S = [T/N|P]-_, % current tally is T/N
        tally_si(T/N),
        (   #T * 6 #> #N,  % current toxicity rate is WORSE THAN 1/6
            length(P, Rec)       % ==> Rec = next-lower dose
        ;   #T * 6 #=< #N, % current toxicity rate is NO WORSE than 1/6
            length([T/N|P], Rec) % ==> Rec = current dose
        )
    ;   S = []-_, % current dose is zero (TODO: does this ever happen?)
        Rec = 0
    ).

%?- stopstate_rec([1/6,0/3]-[2/3], MTD).
%@    MTD = 2
%@ ;  false.

%?- stopstate_rec([T/N,0/3]-[2/3], 2).
%@    T = 0, N = 0
%@ ;  T = 0, N = 1
%@ ;  T = 0, N = 2
%@ ;  T = 0, N = 3
%@ ;  T = 0, N = 4
%@ ;  T = 0, N = 5
%@ ;  T = 0, N = 6
%@ ;  T = 1, N = 6
%@ ;  false.

%% TODO: Note the semantics here are that stopstate(S, MTD) means
%%       that *IF* S is a valid 'stop' state (i.e., one from which
%%       no further dose-escalation decision is permissible),
%%       *THEN* MTD is the recommended dose.
%?- stopstate_rec([T/N,0/3]-[2/3], 1).
%@    T = 1, N = 1
%@ ;  T = 1, N = 2
%@ ;  T = 2, N = 2
%@ ;  T = 1, N = 3
%@ ;  T = 2, N = 3
%@ ;  T = 3, N = 3
%@ ;  T = 1, N = 4
%@ ;  T = 2, N = 4
%@ ;  T = 3, N = 4
%@ ;  T = 4, N = 4
%@ ;  T = 1, N = 5
%@ ;  T = 2, N = 5
%@ ;  T = 3, N = 5
%@ ;  T = 4, N = 5
%@ ;  T = 5, N = 5
%@ ;  T = 2, N = 6
%@ ;  T = 3, N = 6
%@ ;  T = 4, N = 6
%@ ;  T = 5, N = 6
%@ ;  T = 6, N = 6
%@ ;  false.

%% Can I now ensure that any declared MTD obeys the Korn'94 definition?
%% "The MTD is then defined as the highest dose level (≥ 1)
%%  in which 6 patients have been treated with ≤ 1 instance
%%  of DLT, or dose level 0 if there were ≥ 2 instances of DLT
%%  at dose level 1."
%?- setof(StopState, Path^(phrase(path([]-[0/0,0/0]), Path), phrase((..., [StopState,stop,recommend_dose(2)]), Path)), MTD2States).
%@    MTD2States = [[0/6,0/3]-[],[0/6,1/6]-[],[1/6,0/3]-[],[1/6,1/6]-[]].

%?- setof(StopState, Path^(phrase(path([]-[0/0,0/0]), Path), phrase((..., [StopState,stop,recommend_dose(1)]), Path)), MTD1States).
%@    MTD1States = [[0/6]-[2/3],[0/6]-[2/6],[0/6]-[3/3],[0/6]-[3/6],[0/6]-[4/6],[1/6]-[2/3],[1/6]-[2/6],[1/6]-[3/3],[...]-[...],... - ...|...].

%?- setof(StopState, Path^(phrase(path([]-[0/0,0/0]), Path), phrase((..., [StopState,stop,recommend_dose(0)]), Path)), MTD0States).
%@ caught: error(existence_error(procedure,path/3),path/3) % Scryer #1246
%@    MTD0States = [[2/3]-[0/0],[2/6]-[0/0],[2/6]-[2/3],[2/6]-[2/6],[2/6]-[3/3],[2/6]-[3/6],[2/6]-[4/6],[3/3]-[0/0],[...]-[...],... - ...|...].

%% Right away, we can see that 'des' occurs when we should recommend_dose/1.
%% This makes me wonder whether letting 'des' be some kind of default
%% is wrong. Perhaps I should also have a concept of regretting 'des'?
%% Or else maybe there should be some *more active* effort to recommend_dose/1
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
%@    Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[0/6,0/3]-[],stop,recommend_dose(2)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[1/6,0/3]-[],stop,recommend_dose(2)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[0/6]-[2/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[1/6]-[2/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[2/6]-[2/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[2/6,0/3]-[],des,[3/6]-[2/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[0/6]-[3/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[1/6]-[3/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[2/6]-[3/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[0/3,0/3]-[],sta,[3/6,0/3]-[],des,[3/6]-[3/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[1/6,0/3]-[],stop,recommend_dose(2)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[0/6]-[2/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[1/6]-[2/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[2/6]-[2/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[2/6,0/3]-[],des,[3/6]-[2/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[0/6]-[3/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[1/6]-[3/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[2/6]-[3/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[3/6,0/3]-[],des,[3/6]-[3/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[0/6]-[4/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[1/6]-[4/6],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[2/6]-[4/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[1/3,0/3]-[],sta,[4/6,0/3]-[],des,[3/6]-[4/6],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[0/6]-[2/3],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[1/6]-[2/3],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[2/6]-[2/3],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[2/3,0/3]-[],des,[3/6]-[2/3],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[0/6]-[3/3],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[1/6]-[3/3],stop,recommend_dose(1)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[2/6]-[3/3],stop,recommend_dose(0)]
%@ ;  Path = [sta,[0/3]-[0/0],esc,[3/3,0/3]-[],des,[3/6]-[3/3],stop,recommend_dose(0)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[0/6,1/6]-[],stop,recommend_dose(2)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[1/6,1/6]-[],stop,recommend_dose(2)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[2/6,1/6]-[],stop,recommend_dose(1)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[0/3,1/6]-[],sta,[3/6,1/6]-[],stop,recommend_dose(1)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[1/6,1/6]-[],stop,recommend_dose(2)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[2/6,1/6]-[],stop,recommend_dose(1)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[3/6,1/6]-[],stop,recommend_dose(1)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[1/3,1/6]-[],sta,[4/6,1/6]-[],stop,recommend_dose(1)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[2/3,1/6]-[],stop,recommend_dose(1)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[1/6]-[0/0],esc,[3/3,1/6]-[],stop,recommend_dose(1)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[2/6]-[0/0],stop,recommend_dose(0)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[3/6]-[0/0],stop,recommend_dose(0)]
%@ ;  Path = [sta,[1/3]-[0/0],sta,[4/6]-[0/0],stop,recommend_dose(0)]
%@ ;  Path = [sta,[2/3]-[0/0],stop,recommend_dose(0)]
%@ ;  Path = [sta,[3/3]-[0/0],stop,recommend_dose(0)]
%@ ;  false. % J=46 paths!

%% ~~~ BENCHMARKS ~~~

%?- time(J+\(length(Ds,1), maplist(=(0/0), Ds), findall(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J))).
%@    % CPU time: 0.782s % rebis-dev 1/23
%@    J = 10.
%@    % CPU time: 1.140s % rebis-dev
%@    J = 10
%@ ;  % CPU time: 0.000s
%@    false.
%@    % CPU time: 4.148s % master
%@    J = 10.
%?- time(J+\(length(Ds,2), maplist(=(0/0), Ds), findall(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J))).
%@    % CPU time: 3.995s % rebis-dev 1/23
%@    J = 46.
%@    % CPU time: 5.741s % rebis-dev
%@    J = 46
%@ ;  % CPU time: 0.000s
%@    false.
%@    % CPU time: 20.956s % master
%@    J = 46.
%?- time(J+\(length(Ds,3), maplist(=(0/0), Ds), findall(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J))).
%@    % CPU time: 14.185s % rebis-dev 1/23
%@    J = 154.
%@    % CPU time: 20.284s % rebis-dev
%@    J = 154
%@ ;  % CPU time: 0.000s
%@    false.
%@    % CPU time: 75.597s % master
%@    J = 154.
%?- time(J+\(length(Ds,4), maplist(=(0/0), Ds), findall(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J))).
%@    % CPU time: 42.324s % rebis-dev 1/23
%@    J = 442.
%@    % CPU time: 60.256s % rebis-dev
%@    J = 442
%@ ;  % CPU time: 0.000s
%@    false.
%@    % CPU time: 227.109s % master
%@    J = 442.
%?- time(J+\(length(Ds,5), maplist(=(0/0), Ds), findall(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J))).
%@    % CPU time: 115.892s % rebis-dev 1/23
%@    J = 1162.
%@    % CPU time: 162.863s % rebis-dev
%@    J = 1162
%@ ;  % CPU time: 0.000s
%@    false.
%@    % CPU time: 601.521s % master
%@    J = 1162.
%?- time(J+\(length(Ds,6), maplist(=(0/0), Ds), findall(Path, phrase(path([]-Ds), Path), Paths), length(Paths, J))).
%@    % CPU time: 295.461s % rebis-dev 1/23
%@    J = 2890.
%@    % CPU time: 421.905s % rebis-dev
%@    J = 2890
%@ ;  % CPU time: 0.000s
%@    false.
%@    % CPU time: 1529.402s % master
%@    J = 2890.


%% What does a path look like, under optional size-2/3 cohorts?
%?- length(Ds,2), maplist(=(0/0), Ds), phrase(path([]-Ds), Path).
%@    D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[0/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[0/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[0/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[0/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[0/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[1/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[1/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[1/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[1/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[1/4,...]-[],sta,... - ...|...]
%@ ;  D = [0/0,0/0], Path = [esc,[0/2]-[0/0],sta,[0/4]-[0/0],esc,[0/2,0/4]-[],sta,[1/4,...]-[],sta,... - ...|...]
%@ ;  ...

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

