%% Prototyping some checks of CPE ...

%%:- use_module('/Users/david/Precis/precautionary/exec/prolog/tally.pl').
:- use_module('/Users/david/Precis/precautionary/exec/prolog/ccd.pl').

:- use_module(library(lambda)).

%% Read path 'matrices' from a file, to speed testing...
paths_j(Ms, J) :-
    open("testD2", read, Stream),
    read_term(Stream, Ms, []), % Ms were generated with ccd:default_ccd(CCD), D=2.
    length(Ms, J).

%% TODO: Implement checks derived from the condition 𝚺ⱼ 𝜋ⱼ ≡ 1
path_probs_sum_1(Paths, Qjd, Njd, Nj_, MaxN) :-
    %% TODO: Once this matrix notation 'settles down', it might be nice
    %%       to consider how goal expansions might implement a cleaner
    %%       notation.
    % 1. Extract dose-wise tallies as a matrix
    maplist(\Path^Qd^( % (Qd) for a length-D list of tallies Q
		Path = (Ls^Rs^Es~>_), foldl(append, [Ls,Rs,Es], [], Qd)
	    ), Paths, Qjd), % (Qjd) is a row-major J*D matrix (J-list of D-lists)
    % 2. Obtain J*D matrices (N|T|U)dj, dose-wise lists of path-wise (enrollment|tox|no-tox)
    maplist(\Qd^Nd^Td^Ud^(
		maplist(\Q^N^T^U^(
			    Q=T/N, U #= N-T
			), Qd, Nd, Td, Ud) % length-D lists (as if we 'attached' a 'd' index)
	    ), Qjd, Njd, Tjd, Ujd), % now length-J lists of length-D lists
    % 3. Sum over d index for j-indexed (path-wise) totals
    maplist(\Nd^N_^( % does underscore serve well to denote a summed-over index?
		sum(Nd, #=, N_)
	    ), Njd, Nj_), % Nj_ =  𝚺_d Njd ∀ j
    % 4. Obtain maxⱼ(Nⱼ) over all paths (needed below to avoid fractions in test)
    foldl(\A^B^MaxAB^(
	      MaxAB #= max(A,B)
	  ), Nj_, 0, MaxN),
    % 5. Check the basic identity 𝚺ⱼ 2^(MaxN - Nⱼ) = 2^MaxN,
    %    derived by plugging p = q = 1/2 into  𝚺ⱼ 𝛑ⱼ ≡ 1.
    %% Let's call the LHS terms Gdⱼ (G for deGeneracy); then ...
    maplist(MaxN+\N^G^(
		2^(MaxN-N) #= G
	    ), Nj_, Gj), % Gj = 2^(MaxN-Nj) ∀ j
    sum(Gj, #=, 2^MaxN),
    true. /* ~~snip~~
    % 6. Check the identities 0 = 𝚺ⱼ 2^(MaxN - Nⱼ)*(Tⱼd-Uⱼd) ∀ d,
    %    derived by taking derivative 𝜕/𝜕p (𝚺ⱼ 𝛑ⱼ ≡ 1) at p=(1/2,1/2,...1/2).
    */

%?- J^Qjd^Nj_^MaxN+\(paths_j(Ps, J), path_probs_sum_1(Ps, Qjd, Njd, Nj_, MaxN)).
%@    J = 212, Qjd = [[0/6,0/1],[1/6,0/1],[1/6,0/1],[2/6,0/1],[1/6,0/1],[2/6,0/1],[2/6,0/1],[0/2,... / ...],[... / ...|...],...|...], Nj_ = [7,7,7,7,7,7,7,8,12,12|...], MaxN = 12.

%?- t+\(paths_j(Ps, J), time(path_probs_sum_1(Ps, Qjd, Njd, Nj_, MaxN))).
%@    % CPU time: 12.205 seconds
%@    % CPU time: 12.209 seconds
%@    true.
