% Attempt an enumeration of ALL 3+3 trials


% Predefined dose levels
dose(D) :- member(D, [1,2,3,4,5]).

% 0 - 3 DLTs per cohort
tox(N) :- member(N, [0,1,2,3]).

% TODO: Distinguish an MTD-not-found terminal case.

% Trial terminates with declared MTD:
declare_mtd(D, [D ^ 1 | Prev]) :-
    stay(D, Prev).
declare_mtd(D_1, [D ^ N | Prev]) :-
    tox(N), N > 1,
    dose(D),
    D_1 is D - 1,
    escalate(D, Prev).

stay(D, [D ^ 1 | Prev]) :-
    dose(D),
    escalate(D, Prev).

escalate(1, []). % base case for recursion
escalate(D1, [D ^ 0 | Prev]) :-
    dose(D),
    D1 is D + 1,
    dose(D1),
    (   escalate(D, Prev)
    ;   stay(D, Prev)
    ).

