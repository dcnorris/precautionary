%% Benchmarking for ccd module

:- use_module(ccd).

%% Because ccd:regression FAILS to indicate that regression has not occurred,
%% we employ CONJUNCTION (;) in order to run this test as a benchmark, under
%% variously toggled goal expansions.
%%
%% Calls to consult/1 have CUMULATIVE side-effects, despite the conjunction,
%% so only one selected 'elegance' may be restored before restarting session
%% as e.g. with C-u C-u [F10] in ediprolog.
benchmark(Elegance) :-
    (	format("Pure code that will go in the paper:~n", []),
	ccd:regression
    ;	format("Faster code, exploiting several inelegances:~n", []),
	consult(inelegance),
	use_module(ccd), % recompile
	ccd:regression
    ;	format("Partial restoration of elegance ... [~q].~n", [Elegance]),
	call(consult(Elegance)), % restore just this elegance
	consult(inelegance),
	use_module(ccd), % recompile
	ccd:regression
    ).

%?- benchmark(elegant_if_). % remember to use C-u C-u F10
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.794 seconds
%@    % CPU time: 11.800 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.128 seconds
%@    % CPU time: 1.132 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_if_].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 2.260 seconds
%@    % CPU time: 2.264 seconds
%@  J(1) = 20.
%@ false.
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 10.995 seconds
%@    % CPU time: 11.000 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.035 seconds
%@    % CPU time: 1.039 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_if_].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 2.083 seconds
%@    % CPU time: 2.087 seconds
%@  J(1) = 20.
%@ false.

%?- benchmark(elegant_qcompare). % remember to use C-u C-u F10
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.280 seconds
%@    % CPU time: 11.285 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.084 seconds
%@    % CPU time: 1.090 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_qcompare].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 10.436 seconds
%@    % CPU time: 10.440 seconds
%@  J(1) = 20.
%@ false.
%@ Pure code that will go in the paper:
%@  D = 1 ...   % CPU time: 11.039 seconds
%@    % CPU time: 11.044 seconds
%@  J(1) = 20.
%@ Faster code, exploiting several inelegances:
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 1.043 seconds
%@    % CPU time: 1.048 seconds
%@  J(1) = 20.
%@ Partial restoration of elegance ... [elegant_qcompare].
%@ Warning: overwriting goal_expansion/2
%@  D = 1 ...   % CPU time: 9.901 seconds
%@    % CPU time: 9.905 seconds
%@  J(1) = 20.
%@ false.
