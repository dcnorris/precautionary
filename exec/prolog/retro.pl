% prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

/* * * * * * MOTIVATING PRINCIPLE * * * * * *

INTUITION: What if, at each enrollment decision, we actively engaged the
conflict between the THERAPEUTIC vs SCIENTIFIC AIMS of the trial? Perhaps
that is the ULTIMATE CONTENT of a proper DSL in this realm. We have to
define these aims, and how we negotiate the conflict between them.

With each enrollment decision, we ask what dose seems the best THERAPEUTIC
choice FOR THE PATIENT, while at the same time contributing to progress
toward the DRUG-DEVELOPMENT GOAL of the trial.


 - - - - - - OPERATIONAL CONTENT - - - - - -

Demonstrating that the above principle has OPERATIONAL SUBSTANCE requires
showing that the LOGIC OF ESCALATION MAY BE RECOVERED DEDUCTIVELY from a
STATEMENT OF TRIAL GOALS. To this end, I will allow for REPRESENTATION of
a larger set of possible escalation paths, upon which the TRIAL GOALS ACT
AS A SELECTIVE PRINCIPLE.

I want to allow, e.g., an *illogical* 'random-walk' type of escalation to
be GENERATED at some level in the program, so that at a 'higher' level it
may demonstrably be rejected.

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

 - - - - - IMPLEMENTATION STRATEGY(IES) - - - - -

In a similar manner to clpz's "separating modeling from search", I would
like to achieve some separation of concerns that empowers high-level DSL
definition (modeling) of a dose-escalation trial via automatic selective
principles that operate transparently.

Mutual recursion seems an appropriate guidance to this, as foreshadowed
by the alternating functors of my 2-stack 'aliquots' code. In each step
of the trial, we ought to be presented with a 'menu' of all conceivable
choices (some of them bad!) as defined by some DCG. But concomitantly,
these choices are considered in the context of the trial's multi-period
operation, and its ultimate goal. This consideration must (I think) be
inherently recursive, and DYNAMIC PROGRAMMING seems at least a sensible
intutition if not necessarily the manner of implementation.

Taken literally, this DP conception would have MTD-declaration define
the final step of the trial, whereas previous steps are retrospectively
constrained by deductively working backwards toward some fixed starting
condition.

~ A LIKELIHOOD-PRINCIPLE-like constraint ~ 

The ORDER in which the toxicities are observed in any given dose level
ought to be regarded as IRRELEVANT to trial decision-making. That is,
only a TALLY such as 4-1/6 (1 DLT out of 6 enrolled at dose 4) may inform
our evaluation of this dose, or indeed ANY decision-making in the trial.
In addition to its resemblance to the LIKELIHOOD PRINCIPLE (LP), we may
also note this constraint's connection with '(ir)rational dose assignment'
issues as advanced by Wages & Braun (2018) and Wages & Bagely (2019).

~ Separating SAFETY from THERAPEUTIC and SCIENTIFIC AIMS ~

Since SAFETY is so fundamental a principle, we may do well to 'install'
it at the low level of the DCG that generates conceivable escalations.
The more interesting and computationally substantive aspects of dose-
escalation:

 - Is this dose sufficiently THERAPEUTIC?
 - Would this enrollment decision advance us toward an RP2D?

may OTOH be more suited for consideration alongside the indeterminacy
of toxicity outcomes.

 - - - - - - IMPLEMENTATION TACTICS - - - - - -

[THIS PART IS A WORK IN PROGRESS ... like the code below!]

For example, a goal such as:

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
nor conclude with an RP2D. While some of these check just possibly
might work at compile-time, many of them surely will be of the long-
running sort Markus Triska describes.]

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

 * * * * * * * * * * * * * * * * * * * * * * */

% TODO:
% 1/ DCR
% 2. Condense commentary for current relevance
% 3/ Help tallylist_mtd/2 to terminate
% 4/ Deliver good trial sequences from safe_esc//1
% 5. Investigate informational properties at the enrollment margin.

/* - - -


The idea that THE NEXT PATIENT is always treated at what is somehow
thought to be 'A GOOD DOSE' ought to be preserved, however. But how
that 'goodness' is to be conceived deserves some free exploration.
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
   that pursue the aim of 'declaring MTDi'. Some kind of COMMUNICATION
   between these individual-level titrations would be involved.

 - - - - */

%% =======================================================
%% At each dose, our experience is defined by a TALLY.
%% Initially, this is a DLT TALLY giving a simple fraction.

tally(DLTs/Enrolled) :-
    Enrolled in 1..99, %Enrolled in 1..sup,
    0 #=< DLTs, DLTs #=< Enrolled.

%% Tallies at a given dose MAY BE comparable. (We will generally
%% compare an OBSERVED tally with some HYPOTHETICAL tally that serves
%% as a DECISION THRESHOLD. But such distinctions are nowhere manifest
%% in the code.)

:- op(900, xfx, user:(&>=)). % '_ can't be(come) any better than _'
&>=(T1/N1, T2/N2) :-
    tally(T1/N1),
    tally(T2/N2),
    T1 #>= T2 + max(0, N1 - N2).

:- op(900, xfx, user:(&=<)). % '_ can't be(come) any worse than _'
&=<(Q1, Q2) :- Q2 &>= Q1.

%% NB: The INDETERMINACY of toxicity renders some tally pairs INCOMPARABLE:
%?- 1/3 &=< 2/6.
%@ false.
%?- 1/3 &>= 2/6.
%@ false.
%% This is because a 1/3 tally COULD evolve to 1/6 &=< 2/6,
%% but COULD INSTEAD evolve to 3/6 or 4/6, which are &>= 2/6.

%% NB: When tally pairs ARE comparable, however, the conjunction
%%     of &>= with &=< coincides with a STRICT EQUALITY of BOTH
%%     numerator and denominator, and so is well defined:

:- op(900, xfx, user:(&=)).
&=(Q1, Q2) :- Q1 &>= Q2, Q2 &>= Q1.

%% &= as above ==> STRICT EQUALITY of BOTH numerator & denominator:
%?- T1/N1 &= T2/N2, (T1 #\= T2; N1 #\= N2), [T1, T2, N1, N2] ins 0..99, label([T1,T2,N1,N2]).
%@ false.

/* TODO: Consider defining other comparisons, including &\=
         and strict inequalities &< , &> on an AS-NEEDED BASIS.
         Learning where, how and why the various possible tally
         comparisons DO or (especially.) *DO NOT* prove useful
         may itself yield insights. */

%% =======================================================
%% lowtally/1 and toxictally/1 are judgments about tallies
%% which in general the DSL ought to define:

lowtally(T/N) :-
    tally(T/N),
    % Is there some facility that would let me define &=<
    % so that the following reifies straightforwardly?
    % Q &=< 0/3 #<==> B, %% <<== doesn't reify
    % Instead, I am forced to duplicate a little bit of code:
    T #= 0 #/\ N #>= 3 #<==> B, %% LHS is (T/N &=< 0/3)
    lowtally_(B, T/N).

lowtally_(1, _).              % A low tally is &=< 0/3,
lowtally_(0, Q) :- Q &=< 1/6. % or &=< 1/6.

toxictally(Q) :- Q &>= 2/6.

/* NB: Forming such judgments about a dose strictly based
       on the tally AT THAT DOSE is not entirely sound
       from a pharmacological point of view! We may need
       to rectify this conceptual flaw at some point. */

%% =======================================================
%% Whereas individual doses are subject to 'JUDGMENTS' as above
%% (see the above caveat, as well), dose-escalation DECISIONS
%% are made based on TALLY LISTS that fully account for all the
%% assessments made thus far in a trial.

%% I'll let tallies be positional LISTS, effectively defaulting to 0/0:
tallylist_dose_tally(Tallies, D, T/N) :-
    (	nth1(D, Tallies, T/N)
    ;	length(Tallies, L),
	D #> L,
	T/N = 0/0
    ).

/* NB: This 'default' treatment is a bit awkward, and depends
   for its correctness upon escalation NEVER SKIPPING A DOSE. */

%?- tallylist_dose_tally(Tallies, Dose, Tally).
%@ Tallies = [_7434/_7436|_7588],
%@ Dose = 1,
%@ Tally = _7434/_7436 ;
%@ Tallies = [_7586, _7434/_7436|_8998],
%@ Dose = 2,
%@ Tally = _7434/_7436 ;
%@ Tallies = [_7586, _8996, _7434/_7436|_10408],
%@ Dose = 3,
%@ Tally = _7434/_7436 ... ETC.


%% =======================================================
%% Here, the GOAL of dose-finding trials comes into view!
%% This predicate holds when the given TallyList supports
%% the declaration of MTD.

tallylist_mtd(TallyList, MTD) :-
    length(TallyList, _), % Enumerate fairly when TallyList is nonground.
    append(LowTallies, [Q | _], TallyList),
    length(LowTallies, MTD), % NB: enforces MTD >= 0 for nonground MTD
    % TODO: Introduce a safetally/1 predicate to impose
    %       L &=< 1/6 on the final element of LowTallies.
    maplist(lowtally, LowTallies),
    toxictally(Q).

%?- tallylist_mtd([T/N | _], 0), N #=< 6, labeling([up], [N,T]).
%@ T = N, N = 2 ;
%@ T = 2,
%@ N = 3 ;
%@ T = N, N = 3 ;
%@ T = 2,
%@ N = 4 ;
%@ T = 3,
%@ N = 4 ;
%@ T = N, N = 4 ;
%@ T = 2,
%@ N = 5 ;
%@ T = 3,
%@ N = 5 ;
%@ T = 4,
%@ N = 5 ;
%@ T = N, N = 5 ;
%@ T = 2,
%@ N = 6 ;
%@ T = 3,
%@ N = 6 ;
%@ T = 4,
%@ N = 6 ;
%@ T = 5,
%@ N = 6 ;
%@ T = N, N = 6 ;
%% ...
%% RE-CYCLES AS EXPECTED

%?- tallylist_mtd([T/N | _], 0), N #=< 6, labeling([down], [N,T]).
%@ T = N, N = 6 ;
%@ T = 5,
%@ N = 6 ;
%@ T = 4,
%@ N = 6 ;
%@ T = 3,
%@ N = 6 ;
%@ T = 2,
%@ N = 6 ;
%@ T = N, N = 5 ;
%@ T = 4,
%@ N = 5 ;
%@ T = 3,
%@ N = 5 ;
%@ T = 2,
%@ N = 5 ;
%@ T = N, N = 4 ;
%@ T = 3,
%@ N = 4 ;
%@ T = 2,
%@ N = 4 ;
%@ T = N, N = 3 ;
%@ T = 2,
%@ N = 3 ;
%@ T = N, N = 2 ;
%@ T = N, N = 6 ;
%@ T = 5,
%@ N = 6 ....
%% RE-CYCLES AS EXPECTED

/* Importantly, some of the dose-escalations admitted by this relation
   are WRONG ... AND THIS IS GOOD! I want to show that incorporating
   further structure in the problem, by introducing dynamic-programming
   'lookahead' within dose-escalation sequences, excludes these cases. */

/* The following query has the effect of listing all 30 possible toxic
   tallies (here halting escalation at dose 1) with N #=< 9. */
%:- set_prolog_flag(toplevel_print_anon, false).
%@ true.
%?- tallylist_mtd(TallyList, 0), TallyList = [_T/_N | _], _N #=< 9, label([_N,_T]).
%@ TallyList = [2/2] ;
%@ TallyList = [2/3] ;
%@ TallyList = [3/3] ;
%@ TallyList = [2/4] ;
%@ TallyList = [3/4] ;
%@ TallyList = [4/4] ;
%@ TallyList = [2/5] ;
%@ TallyList = [3/5] ;
%@ TallyList = [4/5] ;
%@ TallyList = [5/5] ;
%@ TallyList = [2/6] ;
%@ TallyList = [3/6] ;
%@ TallyList = [4/6] ;
%@ TallyList = [5/6] ;
%@ TallyList = [6/6] ;
%@ TallyList = [3/7] ;
%@ TallyList = [4/7] ;
%@ TallyList = [5/7] ;
%@ TallyList = [6/7] ;
%@ TallyList = [7/7] ;
%@ TallyList = [4/8] ;
%@ TallyList = [5/8] ;
%@ TallyList = [6/8] ;
%@ TallyList = [7/8] ;
%@ TallyList = [8/8] ;
%@ TallyList = [5/9] ;
%@ TallyList = [6/9] ;
%@ TallyList = [7/9] ;
%@ TallyList = [8/9] ;
%@ TallyList = [9/9] ;
%@ TallyList = [2/2, _424] ;
%@ TallyList = [2/3, _424] ...

%% =======================================================
%% Based on the TallyList STATE at any given point in the trial,
%% there will be 0 or more dose levels that 'look safe' to enroll.
%% The safe_nextdose/2 predicate posts a CLP(Z) constraint upon
%% allowable NextDose. This is the key mechanism for creating the
%% options for escalation.

safe_nextdose(TallyList, NextDose) :-
    maxsafe_nextdose(TallyList, MaxNextDose),
    NextDose in 1..MaxNextDose.

maxsafe_nextdose([], 1).
maxsafe_nextdose([Q|Qs], Dose) :-
    reverse([Q|Qs], Rs),
    maxsafe_nextdose_(Rs, Dose).

maxsafe_nextdose_([R|Rs], Dose) :-
    length([R|Rs], Len),
    (	(   R &=< 0/3
	;   R &=< 1/6
	) -> Dose #= Len + 1
    ;	R &=< 1/3 -> Dose #= Len
    ;	maxsafe_nextdose_(Rs, Dose)
    ).

%?- safe_nextdose([0/3,1/6], SafeDose).
%@ SafeDose in 1..3.

%?- safe_nextdose([0/3, 1/3], SafeDose).
%@ SafeDose in 1..2.

%% ===================================================
%% UPDATES to the TallyList STATE are described here.
%% This is the PURELY ARITHMETICAL RELATION between
%% TALLYLISTs and ESCALATION STEPS.

tallylist0_escalation_tallylist(TallyList, [], TallyList).
tallylist0_escalation_tallylist(TallyList0, [Dose - T/N | Cs], TallyList) :-
    length(TallyList0, Dmaxyet), %% Dmaxyet is highest dose tried thus far
    % I have mixed feelings about HARD-CODING A NO-DOSE-SKIPPING CONSTRAINT,
    % but it's hard to imagine that a rigorously LOGICAL approach could ever
    % dispense with this constraint. The fact that it feels so necessary in
    % my logic PROGRAMMING may well be a 'sign'!
    1 #=< Dose, Dose #=< Dmaxyet + 1,
    tally(T/N),
    % TODO: Employ constraint reification to avoid delicate logic here,
    %       or maybe even try reif:if_/3 at long last!
    (	nth1(Dose, TallyList0, T0/N0, Remainder),
	N1 #= N0 + N,
	T1 #= T0 + T,
	nth1(Dose, TallyList_, T1/N1, Remainder)
    ;	append(TallyList0, [T/N], TallyList_),
	% NB: The following goal renders these clauses OPERATIONALLY disjoint
	length(TallyList_, Dose) % w/o dose-skipping, length(TallyList0, Dose-1) holds.
    ),
    tallylist0_escalation_tallylist(TallyList_, Cs, TallyList).

%?- tallylist0_escalation_tallylist([0/3, 1/6], [3 - 0/3], TallyList).
%@ TallyList = [0/3, 1/6, 0/3] ;
%@ false.

%% =======================================================
%% We use this DCG to represent all MERELY SAFE escalations
%% starting from some given TallyList STATE, terminating when
%% an MTD is declared. This DCG will include EXCESSIVELY SAFE
%% escalations that do nothing to 'advance (the scientfic aims
%% of) the trial'.

%%safe_esc(_) --> [freeze]. % halt runaway DCG for debugging
safe_esc(TallyList0) -->
    (	{ tallylist_mtd(TallyList0, MTD) } -> [declare_mtd(MTD)]
    ;	{ safe_nextdose(TallyList0, SafeDose),
	  tally(T/1),
	  labeling([down], [SafeDose,T]),
	  tallylist0_escalation_tallylist(TallyList0,
					  [SafeDose - T/1],
					  TallyList)
	},
	[SafeDose - T/1],
	safe_esc(TallyList) %% TODO: Does this run depth-first?
    ).

%?- phrase(safe_esc([0/3,2/3]), Hmm).
%@ Hmm = [declare_mtd(1)]. %% CORRECT -- this tallylist already supports an MTD.

%?- phrase(safe_esc([0/3,1/3]), Hmm).
%@ Hmm = [2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 2-0/1, 2-0/1, 3-1/1, 3-1/1, declare_mtd(2)] ;
%@ Hmm = [2-0/1, 2-0/1, 2-0/1, 3-1/1, 3-0/1, 3-1/1, declare_mtd(2)] ;
%@ Hmm = [2-0/1, 2-0/1, 2-0/1, 3-1/1, 3-0/1, 3-0/1, 3-1/1, declare_mtd(2)] ;
%@ Hmm = [2-0/1, 2-0/1, 2-0/1, 3-1/1, 3-0/1, 3-0/1, 3-0/1, 3- ... / ..., declare_mtd(...)] ;
%@ Hmm = [2-0/1, 2-0/1, 2-0/1, 3-1/1, 3-0/1, 3-0/1, 3-0/1, 3- ... / ..., ... - ...|...] ...

%?- length(Hmm, 6), phrase(safe_esc([0/3,1/3]), Hmm).
%@ Hmm = [2-0/1, 2-0/1, 2-0/1, 3-1/1, 3-1/1, declare_mtd(2)] ;
%@ Hmm = [2-0/1, 2-0/1, 2-0/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 2-0/1, 1-1/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 2-0/1, 1-1/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 2-0/1, 1-1/1, 1-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 2-0/1, 1-0/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 2-0/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-1/1, 2-1/1, 1-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-1/1, 2-1/1, 1-0/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-1/1, 2-0/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-1/1, 2-0/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-1/1, 2-0/1, 1-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-1/1, 1-0/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-1/1, 1-0/1, 2-1/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-1/1, 1-0/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-1/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-0/1, 2-0/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-0/1, 2-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-0/1, 1-1/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-0/1, 1-1/1, 2-1/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-0/1, 1-1/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [2-0/1, 1-0/1, 1-1/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-0/1, 1-0/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-0/1, 1-0/1, 1-1/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [2-0/1, 1-0/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-1/1, 2-0/1, 2-1/1, 1-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 2-0/1, 2-1/1, 1-0/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [1-1/1, 2-0/1, 2-0/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 2-0/1, 2-0/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 2-0/1, 2-0/1, 1-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 2-0/1, 1-0/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 2-0/1, 1-0/1, 2-1/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [1-1/1, 2-0/1, 1-0/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 2-0/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-1/1, 1-0/1, 2-0/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 1-0/1, 2-0/1, 2-1/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [1-1/1, 1-0/1, 2-0/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-1/1, 1-0/1, 2-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-1/1, 1-0/1, 1-0/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-1/1, 1-0/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 2-0/1, 2-0/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-0/1, 2-0/1, 2-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 2-0/1, 1-1/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-0/1, 2-0/1, 1-1/1, 2-1/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 2-0/1, 1-1/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-0/1, 2-0/1, 1-1/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 2-0/1, 1-0/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 2-0/1, 1-0/1, 1-1/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 2-0/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-1/1, 2-0/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-0/1, 1-1/1, 2-0/1, 2-1/1, 1-0/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-1/1, 2-0/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ Hmm = [1-0/1, 1-1/1, 2-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-1/1, 1-0/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-1/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 2-0/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 2-0/1, 1-1/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 2-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 1-1/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 1-1/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 1-0/1, 2-0/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 1-0/1, 1-1/1, 2-1/1, declare_mtd(1)] ;
%@ Hmm = [1-0/1, 1-0/1, 1-0/1, 1-0/1, 2-1/1, declare_mtd(1)] ;
%@ false.

%?- length(More, _), phrase(safe_esc([0/3,1/3]), [2-0/1, 2-0/1, 2-0/1 | More]).
%@ More = [3-1/1, 3-1/1, declare_mtd(2)] ;
%@ More = [1-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [3-1/1, 3-0/1, 3-1/1, declare_mtd(2)] ;
%@ More = [3-1/1, 2-0/1, 3-1/1, declare_mtd(2)] ;
%@ More = [3-1/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [3-1/1, 1-0/1, 3-1/1, declare_mtd(2)] ;
%@ More = [3-0/1, 3-1/1, 3-1/1, declare_mtd(2)] ;
%@ More = [3-0/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [2-1/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [2-0/1, 3-1/1, 3-1/1, declare_mtd(2)] ;
%@ More = [2-0/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [1-1/1, 3-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [1-1/1, 3-0/1, 1-1/1, declare_mtd(0)] ;
%@ More = [1-1/1, 2-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [1-1/1, 2-0/1, 1-1/1, declare_mtd(0)] ;
%@ More = [1-1/1, 1-0/1, 1-1/1, declare_mtd(0)] ;
%@ More = [1-0/1, 3-1/1, 3-1/1, declare_mtd(2)] ;
%@ More = [1-0/1, 1-1/1, 1-1/1, declare_mtd(0)] ;
%@ More = [3-1/1, 3-0/1, 3-0/1, 3-1/1, declare_mtd(2)] ....


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% FRIENDLY REMINDER that I should concern myself with what
%% advancement in the trial is guaranteed CONDITIONAL ON D,
%% but REGARDLESS OF NEXT T!


%% =======================================================
%% This relation generalizes tallylist_mtd/2 PROSPECTIVELY.
%% True if, starting from TallyList, the further Escalation
%% (2nd arg) yields MTD.
tallylist_escalation_mtd(TallyList, [], MTD) :-
    tallylist_mtd(TallyList, MTD).
tallylist_escalation_mtd(TallyList, [E|Es], MTD) :-
    safe_nextdose(TallyList, NextDose),
    NextDose #>= MTD, % contraction-operator logic!
    E = NextDose - T/1, tally(T/1),
    tallylist0_escalation_tallylist(TallyList, [E], TallyList1),
    tallylist_escalation_mtd(TallyList1, Es, MTD).

%?- length(Esc, 1), tallylist_escalation_mtd([0/3,1/3], Esc, 1).
%@ Esc = [2-1/1] ;
%@ false.

/* NB: tallylist_escalation_mtd/3 appears NOT TO BE NEEDED! */

/*
It seems essential to the intuition of dynamic programming that it
work somehow from the MTD backward.

There is an interesting tension between the principle that decisions
are made path-independently based on tallylists alone, and the practical
requirement for fair enumerations. It might well be that predicates such
as tallylist_mtd/2 cannot by themselves support fair enumeration!

Perhaps I have a complex enough situation that some of these predicates
must eschew the MGQ. It may well be that the dynamic programming itself
aims to provide this.

Alternatively, perhaps I achieve this by constraining the allowable
dose-escalation sequences within the DCG itself. Thus, the DCG becomes
the principle by which fair enumeration is possible.
*/

%?- tallylist_mtd(TL, 1), length(TL,2).
%@ TL = [0/_10828, _10838/_10840],
%@ _10828 in 3..99,
%@ _10838 in 2..99,
%@ _10838#>=_10906,
%@ _10840#>=_10838,
%@ _10906 in 2..95,
%@ 2+_10972#=_10906,
%@ _10972 in 0..93,
%@ _10972#>=_11014,
%@ _10972#=max(0, _11014),
%@ _11014 in -4..93,
%@ _11014+6#=_10840,
%@ _10840 in 2..99 ;
%@ TL = [1/_27954, _27964/_27966],
%@ _27954 in 6..99,
%@ _28012+_27954#=6,
%@ _28012 in -93..0,
%@ _27964 in 2..99,
%@ _27964#>=_28080,
%@ _27966#>=_27964,
%@ _28080 in 2..95,
%@ 2+_28146#=_28080,
%@ _28146 in 0..93,
%@ _28146#>=_28188,
%@ _28146#=max(0, _28188),
%@ _28188 in -4..93,
%@ _28188+6#=_27966,
%@ _27966 in 2..99 ;
%@ Action (h for help) ? abort
%@ % Execution Aborted

%% True if Count is the MINIMUM number of patients who could be
%% enrolled starting from TallyList, to yield declaration of MTD.
%% TODO: Would tabling this enforce a unique minimum?
tallylist_minstepsto_mtd(TallyList, Count, MTD) :-
    (	tallylist_mtd(TallyList, MTD) -> Count #= 0
    ;	tallylist_mtd(TallyList, _) -> false % If any other MTD applies, fail.
    ;	findall(SafeDose, safe_nextdose(TallyList, SafeDose), SafeDoses),
	% At this point, SafeDoses contains all safe doses that we could enroll next.
	% For each one, I want to find the residual minstepsto MTD, and then I want
	% to take the minimum over all these numbers and add 1.
	% Why not employ a DCG at THIS POINT, fairly enumerating on the length?
	% Wouldn't that yield the 'cleanest' expression of the idea?
	true
    ).

%?- tallylist_escalation_mtd(TL, [], 1).
%@ TL = [0/_10710, _10720/_10722],
%@ _10710 in 3..99,
%@ _10720 in 2..99,
%@ _10720#>=_10788,
%@ _10722#>=_10720,
%@ _10788 in 2..95,
%@ 2+_10854#=_10788,
%@ _10854 in 0..93,
%@ _10854#>=_10896,
%@ _10854#=max(0, _10896),
%@ _10896 in -4..93,
%@ _10896+6#=_10722,
%@ _10722 in 2..99 ;
%@ TL = [1/_27836, _27846/_27848],
%@ _27836 in 6..99,
%@ _27894+_27836#=6,
%@ _27894 in -93..0,
%@ _27846 in 2..99,
%@ _27846#>=_27962,
%@ _27848#>=_27846,
%@ _27962 in 2..95,
%@ 2+_28028#=_27962,
%@ _28028 in 0..93,
%@ _28028#>=_28070,
%@ _28028#=max(0, _28070),
%@ _28070 in -4..93,
%@ _28070+6#=_27848,
%@ _27848 in 2..99 ;
%@ TL = [0/_9132, _9142/_9144, _9148],
%@ _9132 in 3..99,
%@ _9142 in 2..99,
%@ _9142#>=_9216,
%@ _9144#>=_9142,
%@ _9216 in 2..95,
%@ 2+_9282#=_9216,
%@ _9282 in 0..93,
%@ _9282#>=_9324,
%@ _9282#=max(0, _9324),
%@ _9324 in -4..93,
%@ _9324+6#=_9144,
%@ _9144 in 2..99 ;
%@ TL = [1/_26270, _26280/_26282, _26286],
%@ _26270 in 6..99,
%@ _26334+_26270#=6,
%@ _26334 in -93..0,
%@ _26280 in 2..99,
%@ _26280#>=_26402,
%@ _26282#>=_26280,
%@ _26402 in 2..95,
%@ 2+_26468#=_26402,
%@ _26468 in 0..93,
%@ _26468#>=_26510,
%@ _26468#=max(0, _26510),
%@ _26510 in -4..93,
%@ _26510+6#=_26282,
%@ _26282 in 2..99 ;
%@ TL = [0/_7842, _7852/_7854, _7858, _7864],
%@ _7842 in 3..99,
%@ _7852 in 2..99,
%@ _7852#>=_7932,
%@ _7854#>=_7852,
%@ _7932 in 2..95,
%@ 2+_7998#=_7932,
%@ _7998 in 0..93,
%@ _7998#>=_8040,
%@ _7998#=max(0, _8040),
%@ _8040 in -4..93,
%@ _8040+6#=_7854,
%@ _7854 in 2..99 ;
%@ TL = [1/_24992, _25002/_25004, _25008, _25014],
%@ _24992 in 6..99,
%@ _25062+_24992#=6,
%@ _25062 in -93..0,
%@ _25002 in 2..99,
%@ _25002#>=_25130,
%@ _25004#>=_25002,
%@ _25130 in 2..95,
%@ 2+_25196#=_25130,
%@ _25196 in 0..93,
%@ _25196#>=_25238,
%@ _25196#=max(0, _25238),
%@ _25238 in -4..93,
%@ _25238+6#=_25004,
%@ _25004 in 2..99 ;
%@ TL = [0/_6900, _6910/_6912, _6916, _6922, _6928],
%@ _6900 in 3..99,
%@ _6910 in 2..99,
%@ _6910#>=_6996,
%@ _6912#>=_6910,
%@ _6996 in 2..95,
%@ 2+_7062#=_6996,
%@ _7062 in 0..93,
%@ _7062#>=_7104,
%@ _7062#=max(0, _7104),
%@ _7104 in -4..93,
%@ _7104+6#=_6912,
%@ _6912 in 2..99 ;
%@ TL = [1/_24062, _24072/_24074, _24078, _24084, _24090],
%@ _24062 in 6..99,
%@ _24138+_24062#=6,
%@ _24138 in -93..0,
%@ _24072 in 2..99,
%@ _24072#>=_24206,
%@ _24074#>=_24072,
%@ _24206 in 2..95,
%@ 2+_24272#=_24206,
%@ _24272 in 0..93,
%@ _24272#>=_24314,
%@ _24272#=max(0, _24314),
%@ _24314 in -4..93,
%@ _24314+6#=_24074,
%@ _24074 in 2..99 ;
%@ TL = [0/_6128, _6138/_6140, _6144, _6150, _6156, _6162],
%@ _6128 in 3..99,
%@ _6138 in 2..99,
%@ _6138#>=_6230,
%@ _6140#>=_6138,
%@ _6230 in 2..95,
%@ 2+_6296#=_6230,
%@ _6296 in 0..93,
%@ _6296#>=_6338,
%@ _6296#=max(0, _6338),
%@ _6338 in -4..93,
%@ _6338+6#=_6140,
%@ _6140 in 2..99 ;
%@ TL = [1/_23302, _23312/_23314, _23318, _23324, _23330, _23336],
%@ _23302 in 6..99,
%@ _23384+_23302#=6,
%@ _23384 in -93..0,
%@ _23312 in 2..99,
%@ _23312#>=_23452,
%@ _23314#>=_23312,
%@ _23452 in 2..95,
%@ 2+_23518#=_23452,
%@ _23518 in 0..93,
%@ _23518#>=_23560,
%@ _23518#=max(0, _23560),
%@ _23560 in -4..93,
%@ _23560+6#=_23314,
%@ _23314 in 2..99 ;
%@ TL = [0/_5466, _5476/_5478, _5482, _5488, _5494, _5500, _5506],
%@ _5466 in 3..99,
%@ _5476 in 2..99,
%@ _5476#>=_5574,
%@ _5478#>=_5476,
%@ _5574 in 2..95,
%@ 2+_5640#=_5574,
%@ _5640 in 0..93,
%@ _5640#>=_5682,
%@ _5640#=max(0, _5682),
%@ _5682 in -4..93,
%@ _5682+6#=_5478,
%@ _5478 in 2..99 ;
%@ TL = [1/_22652, _22662/_22664, _22668, _22674, _22680, _22686, _22692],
%@ _22652 in 6..99,
%@ _22740+_22652#=6,
%@ _22740 in -93..0,
%@ _22662 in 2..99,
%@ _22662#>=_22808,
%@ _22664#>=_22662,
%@ _22808 in 2..95,
%@ 2+_22874#=_22808,
%@ _22874 in 0..93,
%@ _22874#>=_22916,
%@ _22874#=max(0, _22916),
%@ _22916 in -4..93,
%@ _22916+6#=_22664,
%@ _22664 in 2..99 ;
%@ TL = [0/_4932, _4942/_4944, _4948, _4954, _4960, _4966, _4972, _4978],
%@ _4932 in 3..99,
%@ _4942 in 2..99,
%@ _4942#>=_5046,
%@ _4944#>=_4942,
%@ _5046 in 2..95,
%@ 2+_5112#=_5046,
%@ _5112 in 0..93,
%@ _5112#>=_5154,
%@ _5112#=max(0, _5154),
%@ _5154 in -4..93,
%@ _5154+6#=_4944,
%@ _4944 in 2..99 ;
%@ TL = [1/_22130, _22140/_22142, _22146, _22152, _22158, _22164, _22170, _22176],
%@ _22130 in 6..99,
%@ _22224+_22130#=6,
%@ _22224 in -93..0,
%@ _22140 in 2..99,
%@ _22140#>=_22292,
%@ _22142#>=_22140,
%@ _22292 in 2..95,
%@ 2+_22358#=_22292,
%@ _22358 in 0..93,
%@ _22358#>=_22400,
%@ _22358#=max(0, _22400),
%@ _22400 in -4..93,
%@ _22400+6#=_22142,
%@ _22142 in 2..99 ;
%@ TL = [0/_4502, _4512/_4514, _4518, _4524, _4530, _4536, _4542, _4548, _4554],
%@ _4502 in 3..99,
%@ _4512 in 2..99,
%@ _4512#>=_4622,
%@ _4514#>=_4512,
%@ _4622 in 2..95,
%@ 2+_4688#=_4622,
%@ _4688 in 0..93,
%@ _4688#>=_4730,
%@ _4688#=max(0, _4730),
%@ _4730 in -4..93,
%@ _4730+6#=_4514,
%@ _4514 in 2..99 ;
%@ TL = [1/_21712, _21722/_21724, _21728, _21734, _21740, _21746, _21752, _21758, _21764],
%@ _21712 in 6..99,
%@ _21812+_21712#=6,
%@ _21812 in -93..0,
%@ _21722 in 2..99,
%@ _21722#>=_21880,
%@ _21724#>=_21722,
%@ _21880 in 2..95,
%@ 2+_21946#=_21880,
%@ _21946 in 0..93,
%@ _21946#>=_21988,
%@ _21946#=max(0, _21988),
%@ _21988 in -4..93,
%@ _21988+6#=_21724,
%@ _21724 in 2..99 ;
%@ TL = [0/_4140, _4150/_4152, _4156, _4162, _4168, _4174, _4180, _4186, _4192|...],
%@ _4140 in 3..99,
%@ _4150 in 2..99,
%@ _4150#>=_4266,
%@ _4152#>=_4150,
%@ _4266 in 2..95,
%@ 2+_4332#=_4266,
%@ _4332 in 0..93,
%@ _4332#>=_4374,
%@ _4332#=max(0, _4374),
%@ _4374 in -4..93,
%@ _4374+6#=_4152,
%@ _4152 in 2..99 ;
%@ TL = [1/_21362, _21372/_21374, _21378, _21384, _21390, _21396, _21402, _21408, _21414|...],
%@ _21362 in 6..99,
%@ _21468+_21362#=6,
%@ _21468 in -93..0,
%@ _21372 in 2..99,
%@ _21372#>=_21536,
%@ _21374#>=_21372,
%@ _21536 in 2..95,
%@ 2+_21602#=_21536,
%@ _21602 in 0..93,
%@ _21602#>=_21644,
%@ _21602#=max(0, _21644),
%@ _21644 in -4..93,
%@ _21644+6#=_21374,
%@ _21374 in 2..99 ;
%@ TL = [0/_3888, _3898/_3900, _3904, _3910, _3916, _3922, _3928, _3934, _3940|...],
%@ _3888 in 3..99,
%@ _3898 in 2..99,
%@ _3898#>=_4020,
%@ _3900#>=_3898,
%@ _4020 in 2..95,
%@ 2+_4086#=_4020,
%@ _4086 in 0..93,
%@ _4086#>=_4128,
%@ _4086#=max(0, _4128),
%@ _4128 in -4..93,
%@ _4128+6#=_3900,
%@ _3900 in 2..99 ;
%@ TL = [1/_21122, _21132/_21134, _21138, _21144, _21150, _21156, _21162, _21168, _21174|...],
%@ _21122 in 6..99,
%@ _21234+_21122#=6,
%@ _21234 in -93..0,
%@ _21132 in 2..99,
%@ _21132#>=_21302,
%@ _21134#>=_21132,
%@ _21302 in 2..95,
%@ 2+_21368#=_21302,
%@ _21368 in 0..93,
%@ _21368#>=_21410,
%@ _21368#=max(0, _21410),
%@ _21410 in -4..93,
%@ _21410+6#=_21134,
%@ _21134 in 2..99 ;
%@ TL = [0/_3672, _3682/_3684, _3688, _3694, _3700, _3706, _3712, _3718, _3724|...],
%@ _3672 in 3..99,
%@ _3682 in 2..99,
%@ _3682#>=_3810,
%@ _3684#>=_3682,
%@ _3810 in 2..95,
%@ 2+_3876#=_3810,
%@ _3876 in 0..93,
%@ _3876#>=_3918,
%@ _3876#=max(0, _3918),
%@ _3918 in -4..93,
%@ _3918+6#=_3684,
%@ _3684 in 2..99 ;
%@ TL = [1/_20918, _20928/_20930, _20934, _20940, _20946, _20952, _20958, _20964, _20970|...],
%@ _20918 in 6..99,
%@ _21036+_20918#=6,
%@ _21036 in -93..0,
%@ _20928 in 2..99,
%@ _20928#>=_21104,
%@ _20930#>=_20928,
%@ _21104 in 2..95,
%@ 2+_21170#=_21104,
%@ _21170 in 0..93,
%@ _21170#>=_21212,
%@ _21170#=max(0, _21212),
%@ _21212 in -4..93,
%@ _21212+6#=_20930,
%@ _20930 in 2..99 ;
%@ TL = [0/_3488, _3498/_3500, _3504, _3510, _3516, _3522, _3528, _3534, _3540|...],
%@ _3488 in 3..99,
%@ _3498 in 2..99,
%@ _3498#>=_3632,
%@ _3500#>=_3498,
%@ _3632 in 2..95,
%@ 2+_3698#=_3632,
%@ _3698 in 0..93,
%@ _3698#>=_3740,
%@ _3698#=max(0, _3740),
%@ _3740 in -4..93,
%@ _3740+6#=_3500,
%@ _3500 in 2..99 ;
%@ TL = [1/_20746, _20756/_20758, _20762, _20768, _20774, _20780, _20786, _20792, _20798|...],
%@ _20746 in 6..99,
%@ _20870+_20746#=6,
%@ _20870 in -93..0,
%@ _20756 in 2..99,
%@ _20756#>=_20938,
%@ _20758#>=_20756,
%@ _20938 in 2..95,
%@ 2+_21004#=_20938,
%@ _21004 in 0..93,
%@ _21004#>=_21046,
%@ _21004#=max(0, _21046),
%@ _21046 in -4..93,
%@ _21046+6#=_20758,
%@ _20758 in 2..99 ;
%@ TL = [0/_3366, _3376/_3378, _3382, _3388, _3394, _3400, _3406, _3412, _3418|...],
%@ _3366 in 3..99,
%@ _3376 in 2..99,
%@ _3376#>=_3516,
%@ _3378#>=_3376,
%@ _3516 in 2..95,
%@ 2+_3582#=_3516,
%@ _3582 in 0..93,
%@ _3582#>=_3624,
%@ _3582#=max(0, _3624),
%@ _3624 in -4..93,
%@ _3624+6#=_3378,
%@ _3378 in 2..99 ;
%@ TL = [1/_20636, _20646/_20648, _20652, _20658, _20664, _20670, _20676, _20682, _20688|...],
%@ _20636 in 6..99,
%@ _20766+_20636#=6,
%@ _20766 in -93..0,
%@ _20646 in 2..99,
%@ _20646#>=_20834,
%@ _20648#>=_20646,
%@ _20834 in 2..95,
%@ 2+_20900#=_20834,
%@ _20900 in 0..93,
%@ _20900#>=_20942,
%@ _20900#=max(0, _20942),
%@ _20942 in -4..93,
%@ _20942+6#=_20648,
%@ _20648 in 2..99 .
%@ TL = [0/3, 2/2|_8000] ;
%@ TL = [0/3, 2/3|_9196] ;
%@ TL = [0/3, 3/3|_9986] ;
%@ TL = [0/3, 2/4|_11182] ;
%@ TL = [0/3, 3/4|_12004] ;
%@ TL = [0/3, 4/4|_12794] ;
%@ TL = [0/3, 2/5|_13990] ;
%@ TL = [0/3, 3/5|_14812] ;
%@ TL = [0/3, 4/5|_15634] .

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


/*
But clearly, I need some way to generate dose-escalation sequences
while carrying along sufficient context/state to maintain a semblance
of reasonableness.

Getting the balance right, between excessive (random-walk) freedom
in generating these sequences, and premature confinement, will take
some experimentation.

I need to remember that my CORE QUESTION at this stage is how much
can be derived from the DEAL principle. What things *might* I hope
to derive from DEAL? For example, if a dose level already had been
demonstrated to be too toxic, then enrolling additional patients at
that level does not shorten the time to declaring MTD!

I think I want to use a dose_esc DCG only to create a large enough
(but also TOO large!) set of possible future paths for the trial.
It will then be the job of OTHER predicates to refine that large
space of options to a smaller set of 'rational' escalations.

Example heuristics include:

* Always trying the highest dose that would DEFINITELY move the trial
  forward to its conclusion
* Always trying rather the LOWEST such dose
* Always enrolling the dose that lies on the shortest path-to-MTD
* Always enroll the dose lying on the LONGEST path-to-MTD. (In theory,
  this might be justified on grounds that we want to gain sufficient
  experience with the drug at a variety of levels. But it seems to me
  that making this workable would require introducing some additional
  principle that tended to bring the trial to a close.)

If I assume I'm working with SHORTEST-PATH heuristics to begin with,
then I anticipate possibly needing JUST ONE CUT in some predicate.
This is reasonable, since I actually intend to search the solution
space methodically (with fair enumeration, etc.) and indeed define
a correct solution to be the first one! The 'solutions' excluded by
the cut ARE NOT ACTUALLY CORRECT SOLUTIONS!

Once again, note how the SEPARATION OF CONCERNS creates space for
the DSL to define separately ASPECTS OF DOSE-JUDGING and ASPECTS
OF DOSE-ESCALATION-AS-LEARNING.

*/

%% I begin by modeling a dose-escalation ...
%% AHA! It looks to me as if I must carry along the TALLIES as CONTEXT
%% for this DCG. Is this correct? Is there some other way?
%% It seems I am constantly re-learning how to use DCG state vs semicontext!
%% Intuitively -- based on little more than SYNTAX -- I feel that 'lookahead'
%% (maybe with 'replacement') is appropriate for knowing what the last
%% enrolled dose was. OTOH, it does seem as if the current TALLIES have
%% to be managed as STATE.

/*
The FREEDOM preserved here is that this DCG doesn't need to worry about
running away past the point where MTD should have been declared. Halting
is not the responsibility of this DCG. That gets done in another predicate
which I expect will need to employ cut/0.

The RULES expressed in this DCG are much less stringent. For now, I will
say that we don't enroll doses with a toxictally/1, nor with a lowtally/1.

WHAT'S THE CORE INTUITION?

I want this DCG to show all kinds of crazy dose-escalation paths leading
to a given MTD, enumerating them fairly. I should be able to start this
DCG from a 'blank slate' -- i.e., an empty TallyList = []. I should also
be able to start it from 1 enrollment shy of declaring MTD.

What this DCG should do correctly is accumulate the TallyList and keep
the escalation 'connected' (avoid dose-skipping). Fair enumeration is
also important. EVERYTHING ELSE is the responsibility of a supervisory
predicate.

HEY, WHERE'S MY DYNAMIC PROGRAMMING?

I believe it comes in now that tallylist0_escalation_tallylist/3 has been
implemented. By doing a fair enumeration on length of the escalation,
I can stop search at the shortest path to MTD from any given tallylist.

HMMM...

Could I view some of my efforts here as aiming toward FAIR ENUMERATION?
Dose-escalation rules that press 'onward and upward' could be seen as
avoiding unfair enumeration of escalations! Thus, we obtain a purely
logical interpretation of dose-escalation procedures that would normally
be understood from pharmacological or ethical perspectives.

*/
