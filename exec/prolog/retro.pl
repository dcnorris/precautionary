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
    N #<= 6.                    % out of no more than 6 patients.
