% Module of benchmarking utilities

:- module(benchmarking, [
	      version/1,
	      time/2,
	      seconds_since/2,
	      minutes_since/2
	  ]).

version(Version) :-
    '$scryer_prolog_version'(Version).

%% The following adapted from @triska's time:time/1
:- meta_predicate time(0, -).

time(Goal, Elapsed) :-
        '$cpu_now'(T0),
        (   Goal,
            seconds_since(Elapsed, T0)
	;   seconds_since(Elapsed, T0)
	).

since_elapsed(T0, Elapsed) :-
    '$cpu_now'(T),
    Elapsed is T - T0.

seconds_since(Seconds, T0) :-
    '$cpu_now'(T),
    Seconds is T - T0.

minutes_since(Minutes, T0) :-
    seconds_since(Seconds, T0),
    Minutes is Seconds/60.0.
