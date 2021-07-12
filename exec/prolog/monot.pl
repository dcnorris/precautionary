%% Exploration of non-monotonicity in CLP(ℤ)

:- use_module(library(clpz)).

%% I have extracted here the most fundamental relation between tallies,
%% implemented in 'tally.pl' -- namely, REACHABILITY. For example,
%% could a tally of 3/6 (3 toxicities out of 6 patients dosed) ever be
%% ensue during a trial where we've already dosed 4 patients, with zero
%% toxicity? That is the question whether 0/4 ~~ 3/6, which is FALSE.
%% Conversely, 3/6 is a reachable tally from a state of the trial where
%% a 1/4 tally has been observed.
:- op(900, xfx, ~~).
~~(T1/N1, T2/N2) :- % REACHABILITY
    %% The following conditions translate as follows, under 3 scenarios:
    %% (=) If N1==N2, then we get T2 =< T1 =< T2 which holds iff T1==T2.
    %% (<) If N1 < N2, these conds translate to T1 =< T2 =< T1 + MaxTox.
    %% (>) If N1 > N2, these conds translate to T2 =< T1 =< T2 + MaxTox.
    T2 #=< T1 + max(0, N2 - N1),
    T1 #=< T2 + max(0, N1 - N2).


%% As noted above ...

%?- 0/4 ~~ 3/6.
%@ false. % even if patients 5 & 6 have tox, worst tally would be 2/6

%?- 1/4 ~~ 3/6.
%@    true. % if patients 5 & 6 both experience toxicities, we get 1+2=3


%% In investigating reachability, could we ever get tripped up by the
%% non-monotonicity of clpz?

% Could 1/6 ever be reached from an 'earlier' tally T/N with N<6 but T>1?

%?- T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@ false.

% What if I take steps toward indicating my view that N is composed of
% the sum of X's and O's (toxicities and non-toxicities, as commonly
% depicted in dose-escalation diagrams).

%   Added goal
%   vvvvvvvvv
%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@    N = X+O, clpz:(_A#=max(0,_B)), clpz:(_C#=max(0,_D)), clpz:(_B+_E#=6), clpz:(_D+6#=_F), clpz:(T+_A#=_G), clpz:(X+O#=_E), clpz:(X+O#=_F), clpz:(X+O#=_H), clpz:(1+_C#=_I), clpz:(T#=<X+O), clpz:(_A#>=_B), clpz:(_C#>=_D), clpz:(_I#>=T), clpz:(_A in 0..sup), clpz:(_C in 0..sup), clpz:(_I in 2..sup), clpz:(_G in 2..sup), clpz:(T in 2..sup), clpz:(_H in inf..5).

% INTERESTING! I told clpz 'more stuff' (even if I didn't really work
% it out fully, by explaining that O and X are non-negative integers),
% and that left clpz unable to 'tidy up loose ends' during constraint
% propagation?


%?- assertz(clpz:monotonic).
%@    true.

% What if I try that same query again, AFTER making the above assertion?

%?- N = X + O, T/N ~~ 1/6, T #> 1, N #>= T, N #< 6.
%@ caught: error(instantiation_error,instantiation_error(unknown(_1384246),1))

%?- #N = X + O, #T/ #N ~~ 1/6, #T #> 1, #N #>= #T, #N #< 6.
%@ false.

%?- #N = #X + #O, #T/ #N ~~ 1/6, #T #> 1, #N #>= #T, #N #< 6.
%@ false.

%?- #N #= #X + #O, #T/ #N ~~ 1/6, #T #> 1, #N #>= #T, #N #< 6.
%@ false.
