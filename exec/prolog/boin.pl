%% Implement BOIN designs using CCD machinery
:- use_module(library(clpz)).
:- use_module(library(between)).

/*

1. Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials.
   J R Stat Soc C. 2015;64(3):507-523. doi:10.1111/rssc.12089

  Table 1 from [1]:

             cumulative patients treated at current dose (n_j):
              1    2    3    4    5    6    7    8    9    10    11    12

lambda_{1,j} 0/1  0/2  0/3  0/4  0/5  0/6  0/7  1/8  1/9  1/10  1/11  1/12

lambda_{2,j} 1/1  2/2  2/3  2/4  3/5  3/6  4/7  4/8  5/9  5/10  5/11  6/12

elimination   -    -   3/3  3/4  3/5  4/6  4/7  4/8  5/9  5/10  6/11  6/12

*/
% My initial emphasis is on generating all possible paths (CPE) for the BOIN
% design set forth in the table above. Although the BOIN design of [1] lacks
% any terminating principle except elimination of all doses, we do need such
% a rule here. The most natural type of rule, in view of the Table above,
% might be a 'stop-for-consensus' type of rule as found in package 'dtpcrm'.
% This is specified as a maximum number of patients to enroll at any 1 dose.

%% In general, I may require a max-enrollment *parameter*,
%% however unsightly it may be tagging along like this ...
tally(DLTs/Enrolled, MaxN) :-
    ground(MaxN),
    Enrolled in 0..MaxN,
    indomain(Enrolled),
    DLTs in 0..Enrolled,
    indomain(DLTs).

%% Having a 'default' max cohort size will make dev & test a bit easier:
tally(DLTs/Enrolled) :- tally(DLTs/Enrolled, 12).

%% tally/1 terminates:
%?- tally(_), false.
%@ false.

/*

A crucial task for this module is to generate the defining boundaries
for a whole range of BOIN designs. Without access to floating point math,
this involves essentially hard-coding each possible design. An analysis
of the design space thus proves essential.

This also should promote thought about how in general to represent such
design-defining boundaries, and communicate them to the CCD code.

Of course, representing trial designs as Prolog terms (who woulda thot?)
does present itself as one obvious solution. Ideally, this could be done
in such a way that ALL POSSIBLE CCDs could be fairly enumerated.
 
*/

/*

For convenience, we implement the simplest possible BOIN design described
in [1], namely the 'local BOIN' with priors 𝜋_0j = 𝜋_1j = 𝜋_2j, for which
Theorem 2 shows a long-term memory coherence property. Furthermore, we
adopt the recommended values 𝜙1 = 0.6𝜙, and 𝜙2 = 1.4𝜙.

*/

%% These are very close rational approximations to Eqs (2) & (3) in [1],
%% obtained using R MASS function 'fractions':

phipct_lambda1_lambda2(15, 130904/1111271, 5280/29549).

phipct_lambda1_lambda2(20, 723/4598, 224107/939800).

phipct_lambda1_lambda2(25, 9043/45950, 33795/113257).

phipct_lambda1_lambda2(30, 18648/78853, 1172/3269).

%% We also need to tabulate the 5% quantiles of the Beta distribution,
%% from which we will compute the elimination thresholds of [1,p515].

%% Ooh! An interesting question is whether it is logically sound to
%% tabulate a function on an infinite (uncountable!) domain.
%% Perhaps I'm asking, how do I obtain a sound predicate based on
%% the incomplete tabulation of a transcendental function?
%% TODO: Write the elimination-boundary function FIRST, to see what
%%       it tells us it requires in the way of qbeta tabulation.
elim_boundary(T/N, PhiPct) :-
    N in 1..12, indomain(N),
    T in 1..N,
    post05_tally(Qupper, T/N),
    ratless(PhiPct/100, Qupper),
    T_1 #= T - 1,
    post05_tally(Qlower, T_1/N),
    ratless(Qlower, PhiPct/100).

%?- elim_boundary(T/N, 25).
%@    T = 2, N = 2
%@ ;  T = 3, N = 3
%@ ;  T = 3, N = 4
%@ ;  T = 3, N = 5
%@ ;  T = 4, N = 6
%@ ;  T = 4, N = 7
%@ ;  T = 4, N = 8
%@ ;  T = 5, N = 9
%@ ;  T = 5, N = 10
%@ ;  T = 6, N = 11
%@ ;  T = 6, N = 12
%@ ;  false.

%% Less-than relation on rationals
ratless(X1/Y1, X2/Y2) :-
    X1 * Y2 #< X2 * Y1. % 'cross-multiply'

%% The 5% quantile of posterior probability of toxicity,
%% after observing toxicity tally T/N
post05_tally(P, T/N) :-
    Alpha #= T + 1,
    Beta #= N - T + 1,
    qbeta05_alpha_beta(P, Alpha, Beta).

%% qbeta05_alpha_beta(-T/-N, +A, +B) holds if the rational fraction T/N
%% approximates the 0.05 quantile of the Beta(A,B) distribution.
qbeta05_alpha_beta(208010/930249, 2, 1). % 0.223606797749979
qbeta05_alpha_beta(19733/145792, 2, 2). % 0.135350362171584
qbeta05_alpha_beta(21661/58797, 3, 1). % 0.368403149864039
qbeta05_alpha_beta(22399/229471, 2, 3). % 0.0976114628864144
qbeta05_alpha_beta(4365/17558, 3, 2). % 0.248604625730181
qbeta05_alpha_beta(4253/8994, 4, 1). % 0.472870804501588
qbeta05_alpha_beta(2014924/26359415, 2, 4). % 0.0764403914123289
qbeta05_alpha_beta(278748/1472867, 3, 3). % 0.189255377437771
qbeta05_alpha_beta(285655/833806, 4, 2). % 0.342591681998861
qbeta05_alpha_beta(6487/11810, 5, 1). % 0.549280271653059
qbeta05_alpha_beta(5216247/82995322, 2, 5). % 0.0628498917083544
qbeta05_alpha_beta(23734/154961, 3, 4). % 0.153161117975223
qbeta05_alpha_beta(515/1898, 4, 3). % 0.271338372519752
qbeta05_alpha_beta(2574/6155, 5, 2). % 0.418196590747974
qbeta05_alpha_beta(2040/3361, 6, 1). % 0.606962231002917
qbeta05_alpha_beta(827/15494, 2, 6). % 0.0533755004702372
qbeta05_alpha_beta(9139/70979, 3, 5). % 0.128756392804243
qbeta05_alpha_beta(6341/28142, 4, 4). % 0.225321584032448
qbeta05_alpha_beta(3544/10385, 5, 3). % 0.341261436155336
qbeta05_alpha_beta(12409/25890, 6, 2). % 0.479297026408693
qbeta05_alpha_beta(923/1416, 7, 1). % 0.651836344868839
qbeta05_alpha_beta(5669/122205, 2, 7). % 0.0463892639796111
qbeta05_alpha_beta(1/9, 3, 6). % 0.111112706607629
qbeta05_alpha_beta(4254604/22055671, 4, 5). % 0.192902949994131
qbeta05_alpha_beta(8858/30625, 5, 4). % 0.289240816501809
qbeta05_alpha_beta(45623/113969, 6, 3). % 0.400310610809167
qbeta05_alpha_beta(284477/537438, 7, 2). % 0.529320591398868
qbeta05_alpha_beta(6165997/8966688, 8, 1). % 0.687656021933632
qbeta05_alpha_beta(120092/2927419, 2, 8). % 0.0410231675069966
qbeta05_alpha_beta(11687/119564, 3, 7). % 0.0977468134392758
qbeta05_alpha_beta(61667/365433, 4, 6). % 0.168750495987324
qbeta05_alpha_beta(1015256/4038929, 5, 5). % 0.251367627408177
qbeta05_alpha_beta(7736/22427, 6, 4). % 0.344941365943703
qbeta05_alpha_beta(4273/9488, 7, 3). % 0.45035835049339
qbeta05_alpha_beta(45249/79264, 8, 2). % 0.570864452968565
qbeta05_alpha_beta(181177/252733, 9, 1). % 0.716871164436886
qbeta05_alpha_beta(187883/5109482, 2, 9). % 0.0367714378874651
qbeta05_alpha_beta(1662212635/19047996537, 3, 8). % 0.0872644339141503
qbeta05_alpha_beta(107843/718818, 4, 7). % 0.15002824080668
qbeta05_alpha_beta(224371/1008676, 5, 6). % 0.222441101008129
qbeta05_alpha_beta(1742/5739, 6, 5). % 0.303537212564042
qbeta05_alpha_beta(867/2204, 7, 4). % 0.393375783894586
qbeta05_alpha_beta(33653/68248, 8, 3). % 0.493098698936798
qbeta05_alpha_beta(9757/16105, 9, 2). % 0.605836697563495
qbeta05_alpha_beta(52081/70272, 10, 1). % 0.741134449106948
qbeta05_alpha_beta(81278/2439373, 2, 10). % 0.0333192176842298
qbeta05_alpha_beta(149991/1902955, 3, 9). % 0.078820045666029
qbeta05_alpha_beta(252467/1869081, 4, 8). % 0.135075472919641
qbeta05_alpha_beta(57163/286422, 5, 7). % 0.199576149883837
qbeta05_alpha_beta(311744/1149287, 6, 6). % 0.271249914077669
qbeta05_alpha_beta(60509/172976, 7, 5). % 0.349811534571917
qbeta05_alpha_beta(2426/5569, 8, 4). % 0.43562581171077
qbeta05_alpha_beta(54518/102881, 9, 3). % 0.529913200720824
qbeta05_alpha_beta(799/1257, 10, 2). % 0.635640510752758
qbeta05_alpha_beta(18899/24815, 11, 1). % 0.761595809619147
qbeta05_alpha_beta(120602/3959335, 2, 11). % 0.0304601656591372
qbeta05_alpha_beta(1325/18436, 3, 10). % 0.0718702555422698
qbeta05_alpha_beta(59208/481951, 4, 9). % 0.122850663243601
qbeta05_alpha_beta(2406/13291, 5, 8). % 0.181024757242321
qbeta05_alpha_beta(585211/2385697, 6, 7). % 0.245299801274125
qbeta05_alpha_beta(6887/21847, 7, 6). % 0.315237790724452
qbeta05_alpha_beta(9265/23704, 8, 5). % 0.390862301770963
qbeta05_alpha_beta(38443/81331, 9, 4). % 0.472673396439656
qbeta05_alpha_beta(121753/216683, 10, 3). % 0.561894564884689
qbeta05_alpha_beta(282757/427565, 11, 2). % 0.661319331565967
qbeta05_alpha_beta(14649/18803, 12, 1). % 0.779077808054444
qbeta05_alpha_beta(1295/46162, 2, 12). % 0.0280533773201758
qbeta05_alpha_beta(65411/990332, 3, 11). % 0.0660495672163097
qbeta05_alpha_beta(4268416/37885647, 4, 10). % 0.112665780790282
qbeta05_alpha_beta(162358/980071, 5, 9). % 0.165659426715072
qbeta05_alpha_beta(9921/44299, 6, 8). % 0.223955394211319
qbeta05_alpha_beta(134384/468157, 7, 7). % 0.287049002794389
qbeta05_alpha_beta(5515/15544, 8, 6). % 0.354799281019613
qbeta05_alpha_beta(12203/28553, 9, 5). % 0.427380660591787
qbeta05_alpha_beta(13129/25980, 10, 4). % 0.505350269848634
qbeta05_alpha_beta(292/495, 11, 3). % 0.589901395166184
qbeta05_alpha_beta(76961/112572, 12, 2). % 0.68366023523727
qbeta05_alpha_beta(74189513/93416104, 13, 1). % 0.794183334813449
qbeta05_alpha_beta(559411/21516360, 2, 13). % 0.0259993326008678
qbeta05_alpha_beta(6555/107278, 3, 12). % 0.0611029288589776
qbeta05_alpha_beta(50269/483136, 4, 11). % 0.104047307592366
qbeta05_alpha_beta(458/2999, 5, 10). % 0.152717622384663
qbeta05_alpha_beta(15011/72843, 6, 9). % 0.206073335905913
qbeta05_alpha_beta(1397/5300, 7, 8). % 0.26358492378235
qbeta05_alpha_beta(23965/73732, 8, 7). % 0.325028481603775
qbeta05_alpha_beta(2061/5279, 9, 6). % 0.390414867529917
qbeta05_alpha_beta(16891/36720, 10, 5). % 0.459994553476068
qbeta05_alpha_beta(11568/21649, 11, 4). % 0.534343387805782
qbeta05_alpha_beta(2137/3477, 12, 3). % 0.614610317636118
qbeta05_alpha_beta(83920/119329, 13, 2). % 0.703265761051125
qbeta05_alpha_beta(83567/103506, 14, 1). % 0.807363824349865
qbeta05_alpha_beta(7679/316977, 2, 14). % 0.0242257324685366
qbeta05_alpha_beta(1734/30503, 3, 13). % 0.0568468675902468
qbeta05_alpha_beta(46477/480838, 4, 12). % 0.0966583339919568
qbeta05_alpha_beta(178721/1261584, 5, 11). % 0.141663971642153
qbeta05_alpha_beta(10902/57119, 6, 10). % 0.190864686174172
qbeta05_alpha_beta(95576/392143, 7, 9). % 0.243727415766549
qbeta05_alpha_beta(454664/1515615, 8, 8). % 0.29998647413777
qbeta05_alpha_beta(21237/59063, 9, 7). % 0.359565210217654
qbeta05_alpha_beta(118039/279345, 10, 6). % 0.422556337157157
qbeta05_alpha_beta(520423/1063720, 11, 5). % 0.489248110405537
qbeta05_alpha_beta(521/930, 12, 4). % 0.560215564018444
qbeta05_alpha_beta(89213/140149, 13, 3). % 0.636558234458598
qbeta05_alpha_beta(2053/2849, 14, 2). % 0.720603806387123
qbeta05_alpha_beta(1289908/1575049, 15, 1). % 0.818963727477915

%% From these predicates, we should now be able to obtain Prolog terms
%% that define the boundaries of a CCD.
