% Attempt an enumeration of ALL 3+3 trials


% Predefined dose levels
dose(D) :- between(1, 5, D).

% 0 - 3 DLTs per cohort
tox(N) :- between(0, 3, N).

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


% ---- Attempt a more concise program ----

declare(D, [D ^ 1, D ^ 1 | E]) :-
    escalation(D, E).
declare(D_1, [D ^ T | E]) :-
    escalation(D, E),
    tox(T), T > 1,
    D_1 is D - 1.

escalation(1, []).

escalation(D1, [D ^ 0 | E]) :-
    dose(D1),
    D is D1 - 1,
    D > 0,
    escalation(D, E).

% In these goals, should I be unifying E with
% a more detailed representation?

% Do I need 2 different predicates, next/2 and trial/2?
next(D, E) :-
    (	up(D, E)
    ;	sta(D, E)
    ).
up(1, []).
up(U, [D^0 | E]) :- dose(D), U is D+1, next(D,E).
sta(D, [D^1 | E]) :- up(D, E).
mtd(D, [D^1 | E]) :- sta(D, E).
mtd(D_, [D^T | E]) :- tox(T), T>1, up(D, E), D_ is D-1.
result(R, E) :- % TODO: Make #doses a parameter
    (	mtd(D, E), R = mtd(D)
    ;	up(6, E), R = mtd_not_found
    ).

% TODO: Try a formulation via PRODUCTION RULES!
