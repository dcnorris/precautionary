% Attempt an enumeration of ALL 3+3 trials
:- use_module(library(clpfd)).

% Prefix op * for 'generalizing away' goals (https://www.metalevel.at/prolog/debugging)
:- op(920, fy, *). *_.  % *Goal always succeeds

/*

Aim initially to describe the 'Design 1' of Simon &al (1997).
The main difficulty in describing this design seems to be the
arbitrarily distant lookback to earlier cohorts, as required
to ensure that at most 1/6 toxicities occur at the declared
MTD.

What may suffice, however, is to define simple (Markov?) rules,
then PROVE that desired properties (such as '< 1/6') result.

ALTERNATIVELY (and surely *better*!), I might try to state the
desiderata, and allow the dose-escalation process to EMERGE as a
*result* of those specifications. For example, Skolnik &al (2008)
offer this definition:

> The MTD is defined as the dose level at which none or one
> of six participants (0% to 17%) experience a DLT, when at least
> two of three to six participants (33% to 67%) experience a DLT
> at the next highest dose.

To support such a declarative specification, I would have to
implement a scan of the full list generated thus far, perhaps
using semicontext notation as in the tree leaf-counting example
by Triska. The goals at any point would then be:
1. Can I declare an MTD?
2. If not, what should the next cohort be?
3. If I can't create another cohort, then declare "no MTD found".

NOTE HOW WONDERFULLY GOAL-DIRECTED THIS WOULD BE!!

Hmmm!! I think I begin to see implicit concepts REVEALED as I try
to think about the 3+3 in this manner:
a) A dose considered 'too toxic', initially this is (Hi+1)
b) A dose known to be 'safe', initially this is 0

Might I even get away with calling the 'safe' dose MTD?

REFERENCES

Skolnik, Jeffrey M., Jeffrey S. Barrett, Bhuvana Jayaraman, Dimple Patel,
and Peter C. Adamson. “Shortening the Timeline of Pediatric Phase I Trials:
The Rolling Six Design.” JCO 26, no. 2 (Jan 10, 2008): 190–95.
https://doi.org/10.1200/JCO.2007.12.7712.

=====

In all honesty, though, perhaps I should acknowledge that a
'running range' is carried forward implicitly as a state variable,
providing context for (de-)escalation decisions.

esc(D, Lo..Hi) --> ...

*/

sample(D, Set) :-
    D in Set,
    indomain(D).

% Let esc(D, Lo..Hi) *describe* a list of 3+3 cohorts starting from dose D in Lo..Hi
:- discontiguous esc//2.

% Treat escalation past Hi via 'peg' nonterminal.
esc(D, Lo..Hi) --> { D #= Hi + 1 }, peg(Hi, Lo).

%peg(D, _) --> { T in 0..1, indomain(T) }, [D * T], [mtd_notreachedPEG(D)].
%peg(D, Lo) --> { T in 2..3, indomain(T) }, [D * T], { D_1 #= D - 1 }, des(D_1, Lo).

%% NB: The following seems clearer to me than separate goals as above.
%%     But are there any strong reasons to prefer separate goals?
peg(D, Lo) --> [D * T], { T in 0..3, indomain(T) },
	       (   { T in 0..1 } -> [mtd_notreachedPEG(D)]
	       ;   { T in 2..3 } -> { D_1 #= D - 1 }, des(D_1, Lo)
	       ).

% Deal first with case of escalating to the highest allowed dose
esc(D, Lo..Hi) --> [D ^ 0], { D in Lo..Hi, D1 #= D + 1 }, esc(D1, Lo..Hi).
esc(D, Lo..Hi) --> [D ^ 1], { D in Lo..Hi }, sta(D, Lo..Hi).
%% I SUSPECT THE PROBLEM IS IN THIS VICINITY ...
%% Specifically, I see several clauses that allow T in 0..1.
esc(D, Lo..D) --> [D ^ T], { T in 2..3, indomain(T), D_1 #= D - 1 }, des(D_1, Lo).
%% esc(D, Lo..Hi) --> [D ^ T], { D in Lo..Hi, T in 2..3, indomain(T) },
%% 	   (   { T #= 0 } -> (   { D1 #= D + 1 } -> esc(D1, Lo..Hi)
%% 			     ;	 [mtd_notreachedESC(D)]
%% 			     )
%% 	   ;   { T #= 1 } -> sta(D, Lo..Hi)
%% 	   ;   { D_1 #= D - 1 } -> des(D_1, Lo)
%% 	   ).
sta(D, _..D) --> [D - 0], [mtd_notreachedSTA(D)].
sta(D, Lo..D) --> [D - T], { T in 1..3, indomain(T), D_1 #= D - 1 }, des(D_1, Lo).
sta(D, Lo..Hi) --> [D - T], { D in Lo..Hi, D #< Hi, T in 0..3, indomain(T) },
	   (   { T #= 0 } -> (   { D1 #= D + 1, D1 #=< Hi } -> esc(D1, D..Hi)
			     ;	 [declare_mtd(D)]
			     )
	   ;   { D_1 #= D - 1, D_1 #>= Lo } -> [declare_mtd(D_1)]
	   ).
% NB: Implicit in de-escalation to D is that Hi #= D.
des(D, Lo) --> { D #< Lo }, [declare_mtd(D)].
des(D, Lo) --> { D #>= Lo }, [D : T], { T in 0..3, indomain(T) },
	       (   { T in 0..1 } -> [declare_mtd(D)]
	       ;   { D_1 #= D - 1 }, des(D_1, Lo)
	       ).

%% TODO: Write cohorts as (Dose ^ Tox / Enrolled).
%%       This will enable cohorts of size other than 3,
%%       as needed e.g. to describe an accelerated titration phase.
