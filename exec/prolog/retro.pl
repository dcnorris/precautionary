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

 * * * * * * * * * * * * * * * * * * * * * * */

% TODO:
% 1/ DCR
% 2. Condense commentary for current relevance
% 3/ Help tallylist_mtd/2 to terminate
% 4/ Deliver good trial sequences from safe_esc//1
% 5/ Get tallylist_minstepsto_mtd/3 to terminate in general case.
% 6. Investigate informational properties at the enrollment margin.


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
	  labeling([down], [SafeDose]),
	  labeling([up], [T]), % speeds up 'minsteps' search to higher MTDs
	  tallylist0_escalation_tallylist(TallyList0,
					  [SafeDose - T/1],
					  TallyList)
	},
	[SafeDose - T/1],
	safe_esc(TallyList) %% TODO: Does this run depth-first?
    ).

%?- phrase(safe_esc([0/3,2/3]), Hmm).
%@ caught: error(existence_error(procedure,phrase/2),phrase/2)
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
%% Next, we begin a forward-looking set of calculations.

/* I can foresee the likelihood of a whole class of functions
   computable recursively on the tree of possible trial paths,
   capable of supporting decision-making according to Bellman's
   Principle of Optimality.

   In each case, the quantity computed has to be some kind of
   utility or disutility. Two candidates are the maximum and
   minimum duration of the trial, looking forward from any
   given state of the trial. Of these, the MINIMUM clearly
   lends itself more readily to feasible computation, since
   the maximum duration seems to require an exhaustive search
   or else some 'higher' principle along lines of Boolos's
   inference.
*/

%% True if Count is the MINIMUM number of patients who could be
%% enrolled starting from TallyList, to yield declaration of MTD.
%% TODO: Would tabling this enforce a unique minimum, thereby
%%       allowing me to avoid unsound ->? How about reif:if_/3?
%% OTOH: Does this kind of highly selective (extremal) solution-finding
%%       essentially work at the margins of current Prolog capabilities,
%%       or at least at the OUTERMOST LEVEL of the code, and so render
%%       worries about 'soundness' less compelling?
tallylist_minstepsto_mtd(TallyList, MinSteps, MTD) :-
    (	tallylist_mtd(TallyList, MTD) -> MinSteps #= 0
    ;	tallylist_mtd(TallyList, _) -> false % If any other MTD applies, fail.
    ;	length(Escalation, MinSteps), % Starting from MinSteps = 0 and working upward
	% NB: At this point, we actually know that MinSteps=0 won't work.
	%     I should consider working the above 'special cases'
	%     somehow into this general case ..
	MinSteps #>= 1, % ... but pending that I can post this constraint.
	append(Escalation, [declare_mtd(MTD)], Remainder),
	phrase(safe_esc(TallyList), Remainder) -> true
    ).

%% ==========================================================
%% Brief interlude relating to my difficulties using if_ ...

/* This attempted definition didn't work; see below...
tallylist_minstepsto_mtd(TallyList, MinSteps, MTD) :-
    if_(tallylist_mtd(TallyList, MTD),
	MinSteps #= 0,
	(   length(Escalation, MinSteps),
	    append(Escalation, [declare_mtd(MTD)], Fin),
	    phrase(safe_esc(TallyList), Fin)
	)
       ).
*/

%% This section done with Scryer:
%:- use_module(library(reif)).
%@    true.
%?- if_(append([],[],[]), Ans = yes, Ans = no).
%@ caught: error(existence_error(procedure,append/4),append/4)
%@ caught: error(existence_error(procedure,append/3),append/3)
%?- if_(append([],[]), Ans = yes, Ans = no).
%@ caught: error(existence_error(procedure,append/3),append/3)
%% Why does if_/3 believe its first argument is has arity 1 more than its actual?
%?- if_(1 #= 1, Ans = yep, Ans = nope).
%@ ERROR: Unknown procedure: (#=)/3
%@ ERROR:   However, there are definitions for:
%@ ERROR:         (#=)/2
%@ ERROR: 
%@ ERROR: In:
%@ ERROR:   [11] #=(1,1,_61112)
%@ ERROR:   [10] if_(1#=1,_61150=yep,_61156=nope) at /var/folders/pb/d6v8rn4j6x10bzx8qfrs6w6w0000gn/T/ediprolog5788s9J:9
%@ ERROR:    [9] <user>
%@    Exception: (11) #=(1, 1, _61376) ? abort
%@ % Execution Aborted

%% ==================================================
%% Back to tallylist_minstepsto_mtd/3, which works nicely!

%?- tallylist_minstepsto_mtd([0/3,2/3], MinSteps, MTD).
%@ MinSteps = 0,
%@ MTD = 1.

%% In each pair of timings below, we do label([down], SafeDose), and label([up/down], Tox), resp.
%?- time(tallylist_minstepsto_mtd([0/3,1/3], MinSteps, 0)).
%@ % 30,559 inferences, 0.004 CPU in 0.004 seconds (100% CPU, 7183592 Lips) % label([up],[T])
%@ MinSteps = 2.
%@ % 21,027 inferences, 0.002 CPU in 0.002 seconds (99% CPU, 8564969 Lips) % label([down],[T]) ..
%@ MinSteps = 2.
%?- time(tallylist_minstepsto_mtd([0/3,1/3], MinSteps, 1)).
%@ % 5,211 inferences, 0.001 CPU in 0.001 seconds (95% CPU, 5231928 Lips)
%@ MinSteps = 1.
%@ % 3,515 inferences, 0.001 CPU in 0.001 seconds (88% CPU, 4686667 Lips)
%@ MinSteps = 1.
%?- time(tallylist_minstepsto_mtd([0/3,1/3], MinSteps, 2)).
%@ % 312,270 inferences, 0.034 CPU in 0.034 seconds (99% CPU, 9175776 Lips)
%@ MinSteps = 5.
%@ % 279,630 inferences, 0.030 CPU in 0.030 seconds (100% CPU, 9324108 Lips)
%@ MinSteps = 5.
%?- time(tallylist_minstepsto_mtd([0/3,1/3], MinSteps, 3)).
%@ % 12,115,355 inferences, 1.158 CPU in 1.158 seconds (100% CPU, 10464958 Lips)
%@ MinSteps = 8.
%@ % 13,941,261 inferences, 1.330 CPU in 1.331 seconds (100% CPU, 10480693 Lips)
%@ MinSteps = 8.
%?- time(tallylist_minstepsto_mtd([0/3,1/3], MinSteps, 4)).
%@ % 1,227,482,026 inferences, 116.084 CPU in 116.112 seconds (100% CPU, 10574124 Lips)
%@ MinSteps = 11.
%@ % 1,384,418,425 inferences, 130.278 CPU in 130.308 seconds (100% CPU, 10626683 Lips)
%@ MinSteps = 11.
%@ % 1,932 inferences, 0.000 CPU in 0.001 seconds (91% CPU, 3959016 Lips)
%?- time(tallylist_minstepsto_mtd([0/3,1/3], MinSteps, 5)).
%@ Action (h for help) ? abort
%@ % 1,748,576,073 inferences, 183.197 CPU in 185.157 seconds (99% CPU, 9544789 Lips)
%@ % Execution Aborted

%% NB: These 'solutions' are actually WRONG, interpreted straightforwardly.
%%     For example, MTD = 2 is NOT (as one might hope) the value of MTD for
%%     which the minimum number of steps is 8!
%?- tallylist_minstepsto_mtd([0/3,1/3], 8, MTD).
%@ MTD = 2.

%% ==================================
%% TAKING STOCK & LOOKING FORWARD ...


/* - - - -

OBSERVATIONS:

-) We see here at last a NONTRIVIAL calculation being performed.
   It is nontrivial both for its COMPLEXITY and for its CONTENT.
   It would require actual paper-and-pencil combinatorics work
   to figure out some of the answers obtained.

-) I've availed myself of some impure constructs at what SO FAR
   is the outermost level of the code, but won't be for long.
   So it would be nice to use reif:if_/3 where possible.

-) Something that ought to be kept in mind is that the rather strict &__
   inequalities tend to undervalue the information from larger cohorts.
   This tendency 'tips the scales' in favor of changing the dose earlier
   rather than later. I need to be on guard against having hereby imposed
   a tendency that deserved to make a more spontaneous/organic entrance.

   %?- 4/9 &>= 2/6.
   %@ false.

-) This implementation remains FORWARD-LOOKING, and does not exploit
   any OPTIMAL SUBSTRUCTURE of my problem. The timings collected above
   suggest that looking ahead more than 3 doses may be infeasible.
   This points to the need for a forthrightly RETROSPECTIVE programming
   of the solution.

-) OTOH, the (perhaps modest) success of this forward-looking approach
   suggests that some of the basic infrastructure (e.g., tally comparisons)
   will be suitable for a proper application of the BELLMAN PRINCIPLE.

-) The safe_esc//1 DCG seems promising for RETROSPECTIVE work, although
   it doesn't (yet) work in that direction:

   %?- phrase(safe_esc([0/3,1/3]), [C, declare_mtd(1)]).
   %@ C = 2-1/1 ;
   %@ false.

   %?- length(TL,2), phrase(safe_esc(TL), [C, declare_mtd(1)]).
   %@ false.

   %?- tally(Q1), tally(Q2), phrase(safe_esc([Q1,Q2]), [C, declare_mtd(1)]).
   %@ false.

   %% Yet this works fine, provided Q1 and Q2 are instantiated:
   %?- Q1 = 0/3, Q2 = 1/3, tally(Q1), tally(Q2), phrase(safe_esc([Q1,Q2]), [C, declare_mtd(1)]).
   %@ Q1 = 0/3,
   %@ Q2 = 1/3,
   %@ C = 2-1/1 ;
   %@ false.

-) There is another perspective implicit in dose-escalation trials, which
   I have so far neither articulated nor exploited. This is the supposition
   that our LEARNING about dosing MANIFESTS IN THE DOSE ESCALATION ITSELF.
   We could call this DOSE-ESCALATION-AS-LEARNING or 'DEAL'.

-) Seen from the Precautionary Coherence (PC) perspective, of course, this
   is utterly false! But I think that the logic underlying dose-escalation
   IN PRACTICE adheres ~somehow~ to this notion, however bad it may be.
   So REPRESENTING THIS IDEA OBJECTIVELY is useful. I should figure out how
   to achieve this representation, and the present work provides a medium
   for exploring that.

-) The DEAL principle might well provide constraints that render even the
   forward-looking approach adopted here feasible.

-) An early 'win' from lookahead will be showing we can declare MTD upon
   seeing 2/2 toxicities, without filling the full 'cohort of 3' DESPITE
   the fact that the MTD criterion may be EXPRESSED AS 2/3. That is, show
   that the logic of dose-escalation enables tallylist_mtd/2 to 'think
   several moves ahead'.
   This might be as simple as the brute-force step of calculating, from any
   stage of a partially realized escalation, what possible outcomes there
   might be. When the set of possible MTD declarations becomes singular,
   then just declare now.
   That logic may also be extended (probably more feasibly, too) to the
   next enrollment dosing decision. If that decision would be unaffected
   by some given present enrollment decision, then the present decision is
   'wrong' in a sense.

-) One way to conceive of the utility maximized by Bellman Principle is
   as a metric of PROGRESS for a dose-escalation trial, such that the
   next enrollment dose could be chosen on a maximin/minimax criterion
   that let the next dose chosen be the highest one that was sure to
   advance the trial's progress to its SCIENTIFIC GOAL.

-) Example heuristics include:

   * Always trying the highest dose that would DEFINITELY move the trial
     forward to its conclusion
   * Always trying rather the LOWEST such dose
   * Always enrolling the dose that lies on the shortest path-to-MTD
   * Always enroll the dose lying on the LONGEST path-to-MTD. (In theory,
     this might be justified on grounds that we want to gain sufficient
     experience with the drug at a variety of levels. But it seems to me
     that making this workable would require introducing some additional
     principle that tended to bring the trial to a close.)

 - - - - - */
