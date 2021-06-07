% Various goal expansions that offer substantial speedups ...

% Hooks enabling the various expansions to be selectively turned off:
goal_expansion(qcompare_ground_zcompare, true).
goal_expansion(if_lt_via_zcompare, true).

goal_expansion(qcompare(=<, T1/N1, T2/N2, Truth),
	       %% Substitute direct comparisons for reified constraints when possible:
	       (   ground(T1/N1 - T2/N2) ->
		   (   zcompare(C, N2, N1),
		       (   C == (>),
			   (   T1 + N2 - N1 #=< T2 -> Truth = true
			   ;   Truth = false
			   )
		       ;   (C == (=) ; C == (<)),
			   (   T1 #=< T2 -> Truth = true
			   ;   Truth = false
			   )
		       )
		   )
	       ;   % general case (non-ground args 2 or 3) is handled by CLP(Z) ...
		   T1 + max(0, N2 - N1) #=< T2 #<==> B,
		   zo_t(B, Truth)
	       )
	      ) :- qcompare_ground_zcompare.

goal_expansion(if_(X #< X0, Then, Else),
	       (   zcompare(C, X, X0),
		   (   C == (<),
		       Then
		   ;   (C == (=) ; C == (>)),
		       Else
		   )
	       )
	      ) :- if_lt_via_zcompare.
