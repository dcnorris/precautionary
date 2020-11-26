% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

:- use_module(library(clpfd)).

/* - - -

Toward a declarative DSL for dose-escalation trial design...

I'd like to achieve a mutually recursive formulation of dose-escalation,
such that the GOAL of 'declaring an MTD' serves (as much as possible) to
define not only the termination condition, but---recursively---also the
process of enrollment & (de)escalation.

I believe my aim could be described in similar manner to clpz's "separating
modeling from search".

Operationally, my hope is to let MTD-declaration define the final steps
of the trial, while previous steps are 'retrospectively' constrained in
the manner of dynamic programming. For example, a goal such as:

?- mtd_fromescalation_perprotocol(3, Escalation, +Protocol).

ought to return a fairly enumerated set of Escalation terms that describe
the possible paths (sequences of events) that a trial may take if it runs
according to Protocol and declares 3 as the MTD. (The Protocol might be
a term expressed in the DSL. My basic feeling at this point is that a
CONFIGURATION-FILE expressiveness is what we would desire, something akin
to a SAS PROC, say.)

[As an alternative to expressing the Protocol straightforwardly as a
term which the program could INTERPRET, it may be far more desirable
to COMPILE a Prolog program from a DSL Protocol term. This would allow
for various checks to be made, with feedback informing the DSL user
about whether the trial is well-defined (a definite decision arises
in every possible situation) or whether it might somehow get 'stuck'
in a state from which it can neither proceed to enroll another patient
nor conclude with an RP2D.]

This is to some extent an exploration of WHAT KINDS OF GOALS might be
set forth by clinical investigators, such that the types of designs
currently in vogue may be 'recovered' as consequences. Of course, we
are also interested in obtaining HIGHLY ABSTRACT GOALS suited to the
clinical domain and the clinician's---and patient's!---problem.

[Another interesting (but perhaps far-fetched) application might be to
take a given set of decisions (such as mTPI-style decision grids) and
SEARCH FOR a DSL description (or approximation) to the decision grid.
Note that extensions to the DSL that prove necessary to make this work
would themselves be informative even as to the reasonableness of the
decision grids (and underlying models) themselves. There is certainly
a METAGOL spirit to this sort of undertaking.]

What are some basic constraints I am willing to impose from the outset?
I could posit that the ORDER in which the toxicities are observed in any
given dose level is IRRELEVANT to trial decision-making. That is, only
a TALLY such as 4^1/6 (1 DLT out of 6 enrolled at dose 4) informs our
evaluation of this dose, or indeed ANY decision-making in the trial.
(This constraint bears a resemblance to the LIKELIHOOD PRINCIPLE.)

The idea that THE NEXT PATIENT is always treated at what is somehow
thought to be 'A GOOD DOSE' ought to be preserved, however. But how
that goodness is to be conceived deserves some free exploration.
Perhaps there would be some value in allowing (in general) for a whole
range of doses to be considered 'reasonable' at some given time, with
the DSL including language elements supplying a HEURISTIC for choosing
from among a 'reasonable' set containing more than 1 dose.

INVIOLABLE PRINCIPLES:
0. Dose 0 is non-toxic
1. Next patient always treated at a 'therapeutic' dose. This is a
   CONSTRAINT on the dose-escalation PROCESS.
2. The GOAL of the trial is to 'declare' an RP2D of some kind.
3. We adopt at least a HEURISTIC (or search strategy) that somehow
   aims the dose-escalation toward higher doses AS PERMITTED.
4. Toxicities must (in the MATURE DSL) be addressed as ORDINAL 0..5,
   in order to require the user to deal with the real possibilities
   of fatalities (5) and severe toxicities (4) on trial. (While DEVELOPING
   the DSL, however, I will permit myself to address toxicities as binary.)

FREEDOMS TO PRESERVE:
1. Cohort size should be a free parameter. My approach to preserving
   this freedom will be to design the DSL initially around 1-at-a-time,
   sequential {enrollment, dosing, assessment & decision-making}. Any
   larger cohort sizes should either appear as explicit constraints
   imposed through the DSL, or else (preferably!) as CONSEQUENCES
   OF MORE GENERAL CONSIDERATIONS -- e.g., the activation of decision
   points only at certain denominators which would give rise to
   'cohorts' on a 'WLOG'-type principle.
2. Ordinal toxicities (and their underlying logic) demand support
   within this DSL. The language *cannot* be tied to binary 'DLTs',
   since even the most basic--and common!--generalization in practice
   (accelerated titration) depends upon ordinal toxicities. Moreover,
   even in designs ostensibly built upon binary DLTs, there remains
   an implicit distinction in the 'exceptional' treatment given to
   Grade 5 (fatal) toxicities. Indeed, once past the first prototype,
   it may well be desirable to FORCE the user to speak in terms of
   CTCAE toxicity grades 0..5, so that any dichotomization (DLT=T/F)
   becomes conspicuous in the DSL code.
3. At the very least the POSSIBILITY of allowing for dose TITRATION
   should be entertained. This might be handled most elegantly by
   treating a 'dose-titration trial' as an ensemble of N-of-1 trials
   that pursue the aim of 'declaring MTDi'. Some kind of communication
   between these individual-level titrations would be involved.

 - - - - */

/*
I have some hope for recovering the LOGIC OF ESCALATION *deductively*
from A STATEMENT OF TRIAL GOALS. To this end, I will allow for a larger
set of possible escalation paths in my representation. I want to allow,
e.g., an *illogical* 'random-walk' type of escalation to be REPRESENTED
so that the program logic can demonstrate its ability to exclude such
paths from the solution.

An example of 3+3-PROTOCOL-VIOLATING escalation might be as follows...

[1^0/4, 2^0/2, 1^1/2, 2^0/1, 3^2/3]

(The story might be that an extra 4th patient got enrolled at dose1 by
administrative mix-up, and this was realized after 2 patients had already
been enrolled at dose2. So enrollment at dose2 was halted in order to
correct this error by enrolling a 'round figure' of net 6 patients.
After that, the trial returned to 'filling up' a first cohort of 3 patients
at dose2, and then proceeding 'per protocol' (after 2^0/3) to enroll dose3.)

I offer that story only to exhibit a 'protocol violation' that we would
like to RECOGNIZE as a violation by suitable definition & understanding
of the trial's GOALS and CONSTRAINTS.

INTUITION: What if, at each enrollment decision, we actively engaged the
conflict between the THERAPEUTIC vs SCIENTIFIC AIMS of the trial? Perhaps
that is the ULTIMATE CONTENT of a proper DSL in this realm. We have to
define these aims, and how we negotiate the conflict between them.

With each enrollment decision, we ask what dose seems the best THERAPEUTIC
choice FOR THE PATIENT, while at the same time contributing to progress
toward the DRUG-DEVELOPMENT GOAL of the trial.

*/

/*
Here, I'm representing Escalation by a list of Dose-Assessment pairs
for INDIVIDUAL PATIENTS.
For the long term, this enables an ordinal assessment to be employed,
although initially this is just the type Assessment := OKAY | DLT.
*/
dose_escalation(Dose, Escalation) :-
    onpath_toMTD([Dose - okay | Escalation], _), % Enrollment should be part of a
    onpath_toMTD([Dose - dlt | Escalation], _),  % conceivable 'path to RP2D'.
    gooddose_incontext(Dose, Escalation). % But also should look 'therapeutic'.

onpath_toMTD(Escalation, MTD) :-
    (	declare_mtd(Escalation, MTD)      % Either we're already there,
    ;	dose_escalation(E, Escalation),   % or else a further step
	onpath_toMTD([E|Escalation], MTD) % gets us /closer/.
    ). %% TODO: Articulate 'closer' by a METRIC, so escalation becomes a CONTRACTION OPERATOR

%% Dose looks therapeutic (NB: I'm being careful to avoid a word like 'best'
%% which would imply a claim of optimality) after having seen Escalation.
%% ** DETAILS HARD-CODED HERE SHOULD BE EXPRESSED VIA DSL **
%% NB: We want the constraint N>6 to be IMPLICIT, emerging as a CONSEQUENCE.
gooddose_incontext(Dose, Escalation) :-
    tally_fordose_inescalation(T/N, Dose, Escalation),
    N #=< 6, T #< 1,
    D_1 #= Dose - 1, %% Adjacent lower dose showed...
    tally_fordose_inescalation(T_1/N_1, D_1, Escalation),
    (	T_1 #= 0, N_1 #>= 3 % 0/3 or else
    ;	T_1 #= 1, N_1 #>= 6 % 1/6 toxicities.
    ),
    D1 #= Dose + 1, %% Next-higher dose is either...
    tally_fordose_inescalation(T1/N1, D1, Escalation),
    (	N1 #= 0           % thus far untried,
    ;	N1 #> 0, T1 #>= 2 % or else looks too toxic.
    ).

tally_fordose_inescalation(0/sup, 0, _). %% A zero dose is 'supremely nontoxic'
tally_fordose_inescalation(T/N, Dose, Escalation) :-
    tally_fordose_inescalation_acc(T/N, Dose, Escalation, 0/0).

tally_fordose_inescalation_acc(T/N, _, [], T/N).
tally_fordose_inescalation_acc(T/N, Dose, [Dose^T1/N1 | Escalation], T0/N0) :-
    Tacc #= T0 + T1,
    Nacc #= N0 + N1,
    tally_fordose_inescalation_acc(T/N, Dose, Escalation, Tacc/Nacc).
tally_fordose_inescalation_acc(T/N, Dose, [Dose1^_ | Escalation], T0/N0) :-
    Dose1 #\= Dose,
    tally_fordose_inescalation_acc(T/N, Dose, Escalation, T0/N0).

%% I'll treat tallies as lists of terms of the form Tox/Enrolled,
%% in which position tells the dose level.

mtd_escalation(MTD, Escalation) :-
    escalation_tallies(Escalation, Tallies),
    (	nth1(MTD, Tallies, 0/6) % The MTD is a dose we've tried in
    ;	nth1(MTD, Tallies, 1/6) % at least 6 patients, with at most 1 DLT..
    ),
    nth1(MTD+1, Tallies, T/_),  % ..AND for which the next-higher dose
    T #> 2,                     % has shown at least 2 toxicities
    N #=< 6.                    % out of no more than 6 patients.

/* - - - -

Moving past pseudocode, let's try to implement these ideas in a 2-dose setting.

1. First thing I could use is a DCG to generate all possible enrollment strings.

This ought to represent only the most elemental aspects of dose-escalation.
For example, the fact that the next dose is always at most 1 higher than
the present dose -- or, to begin this exploration, that we always move up
OR down by 1 dose.

2. Superimposed on this may be some sense that our LEARNING (i.e., the trial's
PROGRESS toward its DRUG-DEVELOPMENT AIMS) manifests solely in the up-and-down
of dose-escalation itself.

Seen from a PC perspective, of course, this is utterly false! But I do believe
that the logic underlying dose-escalation IN PRACTICE adheres to this idea,
however bad that may be. So REPRESENTING THIS IDEA OBJECTIVELY is useful.

===
Here's the SCIENTIFIC QUESTION I see motivating the next step of this effort:
How far toward a full design can we proceed, simply by taking a few notional
examples of dose-escalation decisions, and 'processing' them through the logic
of DOSE-ESCALATION-AS-LEARNING?

I'm actually going to call that 'DEAL' ;^), and perhaps it turns out to be a
valuable contribution of this effort that it explores the consequences and
limitations of the DEAL principle.

 - - - - -*/

%% I begin by modeling a dose-escalation ...
%% AHA! It looks to me as if I must carry along the TALLIES as CONTEXT
%% for this DCG. Is this correct? Is there some other way?
%% It seems I am constantly re-learning how to use DCG state vs semicontext!
%% Intuitively--based on the SYNTAX, really!--, I feel that 'lookahead'
%% (maybe with 'replacement') is appropriate for knowing what the last
%% enrolled dose was. OTOH, it does seem as if the current TALLIES have
%% to be managed as STATE.
dose_esc([0/0, 0/0, 0/0]) --> [1 - Tox], { Tox in 0..1 }.
dose_esc(Tallies_1), [Dose_1 - Tox_1] -->
                     [Dose_1 - Tox_1, Dose - Tox],
		     { length(D, Tallies_1),
		       Dose in 1..D,
		       Dose #>= Dose_1 - 1, % Next dose in DE always stays
		       Dose #=< Dose_1 + 1, % within 1 level of previous dose.
		       Tox in 0..1 }, % TODO: Eventually, CTCAE Grades 0..5
		     tallies0_assessment_tallies(Tallies_1, Dose - Tox, Tallies),
		     dose_esc(Tallies).
dose_esc(Tallies) --> [declare_mtd(MTD)],
		      { mtd_tallies(MTD, Tallies) }.

%% 'Union' or 'sum' relation on tallies lists: Tallies0 + Tox@Dose = Tallies.
/*
This has become a problem almost immediately! I ought to be able to access
and set elements of a tallies list without so much trouble!
TODO: Switch to something like keypairs lists?
 */
tallies0_assessment_tallies(Tallies0, Dose - Tox, Tallies) :-
    (	nth1(Dose, Tallies0, T0/N0, Remainder),
	N #= N0 + 1,
	T #= T0 + Tox,
	nth1(Dose, Tallies, T/N, Remainder)
    ;	append(Tallies0, [Tox/1], Tallies),
	length(Tallies, Dose) % w/o dose-skipping, length(Tallies0, Dose-1) holds.
    ).

%?- tallies0_assessment_tallies([0/0,0/0,0/0], 2 - 0, Tallies).
%@ Tallies = [0/0, 0/1, 0/0] ;
%@ false.
%@ Tallies = [0/0, 1/1, 0/0] ;
%@ false.
%@ Tallies = [_6368, 1/1|_6376] ;
%@ Action (h for help) ? Unknown option (h for help)
%@ Action (h for help) ? abort
%@ % Execution Aborted
%@ ERROR: Unknown procedure: tallies0_assessment_tallies/3 (DWIM could not correct goal)

%% I'll let tallies be positional LISTS, effectively defaulting to 0/0:
tallies_dose_tally(Tallies, D, T/N) :-
    (	nth1(D, Tallies, T/N)
    ;	length(Tallies, L),
	D #> L,
	T/N = 0/0
    ).

%% MGQ
%?- tallies_dose_tally(Tallies, Dose, Tally).
%@ Tallies = [_2884/_2886|_2898],
%@ Dose = 1,
%@ Tally = _2884/_2886 ;
%@ Tallies = [_2896, _2884/_2886|_4332],
%@ Dose = 2,
%@ Tally = _2884/_2886 ;
%@ Tallies = [_2896, _4330, _2884/_2886|_5766],
%@ Dose = 3,
%@ Tally = _2884/_2886 ;
%@ Tallies = [_2896, _4330, _5764, _2884/_2886|_7200],
%@ Dose = 4,
%@ Tally = _2884/_2886 ... YEP!

/*
TODO: Show that it's possible to declare MTD upon seeing 2/2 toxicities,
      without filling the full 'cohort of 3' DESPITE the fact that the
      MTD criterion may be EXPRESSED AS 2/3. That is, show that the logic
      of dose-escalation enables the DSL to 'think several moves ahead'.

This might be as simple as the brute-force step of calculating, from any
stage of a partially realized escalation, what possible outcomes there
might be. When the set of possible MTD declarations becomes singular,
then just declare now.

That logic may also be extended (probably more feasibly, too) to the next
enrollment dosing decision. If that decision would be unaffected by some
given present enrollment decision, then the present decision is 'wrong'
in a sense.

*/

/* - - - - -

AHA! Here is the RECURSIVE concept I have been waiting for! I ought to
calculate the number of possible escalations FORWARD FROM A GIVEN TALLY
that results in a given outcome, such as declare_mtd(3).

This quantity can be calculated RECURSIVELY and (probably) allows injection
of DYNAMIC-PROGRAMMING ideas. With any luck (or maybe a lot of it!), the
dynamic-programming approach (plus tabling?) may allow me to circumvent an
otherwise impure, forall-based approach to this problem.

What does this circumvention depend on? Ultimately, it depends on a base
case that can be aggregated through the search tree.

PLAN:
/1 Implement stricter tally comparisons &<, &>, &=, etc.
   where the comparison is conceptually grounded in the
   sequences of (0/1 | 1/1) that sum up to them.
2) Implement relation {TallyList - MTD}
3) Relate TallyLists to Escalations that generate them.
   Initially, this can be *ANY* Escalation sequence consistent
   with some DCG. I can leave 'for later' the culling of these
   many candidate Escalations, according to DEAL and other
   principles.
4) Implement relation {TallyList - Escalation - MTD}
5) Render an 'abstract interpretation' of the above relation,
   concerned strictly with the LENGTH of the Escalation needed
   to reach the given MTD, or at least SOME QUANTITY THAT CAN
   BE CALCULATED RECURSIVELY on a tree. This could be, e.g.,
   the minimum number of patients needed to enroll to reach
   an MTD declaration.

I have to remember the IMPORTANCE of such a recursive calculation
derives from its potential to accomplish a forall/3 or \+ BY PURE
MEANS!

What happens at the point of IMPLEMENTATION is the {TallyList - MTD}
relation serves as the base case for {TallyList - Escalation - MTD},
where Escalation = []. Another case deals with the VERY FINITE case
where Escalation = [Dose - Tox] for a single patient. I can probably
deal with this 'in closed form' (in some sense I don't yet understand
fully) so that I can obtain a number. At each node running backwards
toward the root of the tree, I will obtain a new summary. The forall
or \+ happens when I find I've counted ZERO THINGS HAPPENING, say.

It might even be the case that I will require MULTIPLE BINARY OPS
in order to fully characterize the 'answer to the question' at any
given point during the trial. Maybe {max, min, and} are all needed!
Ideally, I would parametrize such ops, for accumulation within the
tree, by a generic maplist-like predicate.

** THE MOST IMPORTANT THINKING to do at this point is to understand
** what sort of \+iness I need to implement in 'deducing' the trial.

CAPTURING *PROGRESS* VIA DYNAMIC PROGRAMMING

I think that DP promises to make a contribution mainly in defining
a notion of 'progress' for dose-escalation trials, such that next
enrollment dose could be chosen on a maximin or minimax criterion
that let the next dose chosen be the highest one that was sure to
reduce the number of steps to MTD.

 - - - - - - */

% This is true if, starting from Tallies a *further* Escalation sequence
% starting from Dose results in declaration of MTD.

% The base case for this recursion defines what tallies yield an MTD right away.
% Note that Dose is anonymous in this clause because MTD declaration depends
% only on the net Tallies seen so far.
/*
tallies_dose_escalation_mtd(Tallies, _, [], MTD) :-
    true.
tallies_dose_escalation_mtd(Tallies, Dose, Escalation, MTD) :-
    true.
*/

/*
Something that ought to be kept in mind is that these rather strict
inequalities tend to undervalue the information from larger cohorts.
This tendency 'tips the scales' in favor of changing the dose earlier
rather than later. I need to be on guard against having hereby imposed
a tendency that deserved to make a more spontaneous/organic entrance.
*/
%?- 4/9 &>= 2/6.
%@ false.

%?- maplist(lowtally, []).
%@ true.

%?- lowtally(0/3).
%@ false.

%?- tally(0/3).
%@ true.

%?- 0/3 &=< 0/3.
%@ true ;
%@ false.

%?- 0/3 &=< 0/3.
%@ true ;
%@ false.



%?- lowtally_(1, 0/3).
%@ true.
%@ false.
%@ true.

%?- maplist(lowtally, [0/3, 0/6]).
%@ false.
%@ false.
%@ true ;
%@ false.
%@ true ;
%@ true.

%?- T/N &>= 1/3.
%@ T in 1..sup,
%@ N#>=T,
%@ N in 1..sup.

%?- 1/10 &>= 1/3.
%@ false.

%?- tallylist_mtd(Tallies, 0).
%@ Tallies = [_928590/_928592|_928586],
%@ _928590 in 1..3,
%@ _928592#>=_928590,
%@ _928592 in 1..3 ;
%@ Tallies = [_931550/_931552|_931546],
%@ _931550 in 2..sup,
%@ _931552#=<_931550+2,
%@ _931552#>=_931550,
%@ _931552 in 4..sup ;
%@ Tallies = [_936246/_936248|_936242],
%@ _936246 in 2..6,
%@ _936248#>=_936246,
%@ _936248 in 2..6 ;
%@ Tallies = [_939206/_939208|_939202],
%@ _939206 in 3..sup,
%@ _939208#=<_939206+4,
%@ _939208#>=_939206,
%@ _939208 in 7..sup ;
%@ Action (h for help) ? abort
%@ % Execution Aborted

%?- 0/3 &>= 0/99.
%@ true.

%?- 0/99 &=< 0/3.
%@ true.

%?- append(L, M, [1,2,3,4,5]).
%@ L = [],
%@ M = [1, 2, 3, 4, 5] ;
%@ L = [1],
%@ M = [2, 3, 4, 5] ;
%@ L = [1, 2],
%@ M = [3, 4, 5] ;
%@ L = [1, 2, 3],
%@ M = [4, 5] ;
%@ L = [1, 2, 3, 4],
%@ M = [5] ;
%@ L = [1, 2, 3, 4, 5],
%@ M = [] ;
%@ false.

%?- maplist(#<(1), [2,3,4]).
%@ true.
%@ false.
%@ ERROR: Syntax error: Unbalanced operator
%@ ERROR: maplist(1 #
%@ ERROR: ** here **
%@ ERROR: <, [2,3,4]) . 

% What does a DLT tally look like?
tally(DLTs/Enrolled) :-
    Enrolled in 1..sup,
    0 #=< DLTs, DLTs #=< Enrolled.

%?- tally(T/N).
%@ T in 0..sup,
%@ N#>=T,
%@ N in 1..sup.

%% Ooh.. look at this! I need to delve a bit into the details
%% of comparing tallies. I rather think that some comparisons
%% ought to be recognized as indeterminate. Perhaps a tally
%% such as 0/2 is INCOMMENSURATE with 1/3! How do you know
%% that the 0/2 won't become a 1/3 with the next patient?
/*
Would the SETS of dose-escalations that correspond to these
tallies provide a sounder basis for thinking about such
comparisons?

Or maybe worst-case scenarios? I can say that 1/2 &=< 2/3
because the possible outcomes are 1/2 --> {2/3, 1/3}, a set
that is at worst equal to 2/3 and possibly better.

Do the complementary-looking comparisons such as [&< and &>=]
lose some of their complementarity under such an understanding?
Or do I merely eject many pairs of tallies from all comparisons?
*/
%?- 0/2 &< 1/3.
%@ true.

:- op(900, xfx, user:(&>=)). % '_ can't be(come) any better than _'
&>=(T1/N1, T2/N2) :-
    tally(T1/N1),
    tally(T2/N2),
    DT #= max(0, N1 - N2),
    T1 #>= T2 + DT.

:- op(900, xfx, user:(&=<)). % '_ can't be(come) any worse than _'
&=<(Q1, Q2) :- Q2 &>= Q1.

%?- 0/2 &>= 0/3.
%@ true.

%?- X #= max(0,2).
%@ X = 2.

%?- T/N &>= 1/3, T in 0..5, N in 0..5.
%@ T in 1..5,
%@ N in 0..3 ;
%@ T in 2..5,
%@ N#=<T+2,
%@ N in 4..5.

%% Define &= as conjunction of &>= with &=<.
:- op(900, xfx, user:(&=)).
&=(Q1, Q2) :- Q1 &>= Q2, Q2 &>= Q1.

%% Demonstrate that &= as above requires STRICT EQUALITY
%% of BOTH numerator & denominator:
%?- T1/N1 &= T2/N2, (T1 #\= T2; N1 #\= N2), [T1, T2, N1, N2] ins 0..6, label([T1,T2,N1,N2]).
%@ false.

%% TODO: Consider defining other comparisons, including &\=
%%       and strict inequalities &< and &>. Learning where,
%%       how and why the various possible tally comparisons
%%       do or--especially!--DO NOT prove useful may itself
%%       yield insights.

lowtally(Q) :-
    tally(Q),
    %(	Q &=< 0/3
    %;	Q &=< 1/6
    %).
    Q &=< 0/3 #<==> B,
    lowtally_(B, Q).

lowtally_(1, _).
lowtally_(0, Q) :- Q &=< 1/6.

tallylist_mtd(Tallies, MTD) :-
    append(LowTallies, [T/N | _], Tallies),
    length(LowTallies, MTD),
    maplist(lowtally, LowTallies),
    tally(T/N),
    (	T/N &>= 1/3
    ;	T/N &>= 2/6
    ).

