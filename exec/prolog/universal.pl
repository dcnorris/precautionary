%% Explore universality of CCD
:- use_module(library(clpz)).
:- use_module(library(lambda)).
:- use_module(library(format)).

:- use_module(ccd).

/*

Initially, let's examine the apparently obvious connection between
the traditional 3+3 design and CCD. If I am not wrong, a CCD with
cohort size 3, cohort limit 6, and the following boundaries yields
our usual 3+3:

             cumulative patients treated at current dose (n_j):
              1    2    3    4    5    6 

lambda_{1,j}  -    -   0/3  0/4  0/5  1/6

lambda_{2,j}  -   2/2  2/3  2/4  2/5  2/6

elimination   -   2/2  2/3  2/4  2/5  2/6


Note that the elimination and de-escalation boundaries are identical!

*/

traditional33(ccd([2/6], [2/6], [0/3, 1/6], 3, 6, 9999)). % 9999 ≈ sup

traditional_d_matrix(D, Matrix) :-
    traditional33(T33),
    ccd_d_matrix(T33, D, Matrix).

%?- traditional33(T33).
%@    T33 = ccd([2/6],[2/6],[0/3,1/6],3,6,9999).

%?- Matrix+\(traditional_d_matrix(2, Matrix)).
%@    Matrix = ([0/3]^[0/6]^[]~>2)
%@ ;  Matrix = ([0/3]^[1/6]^[]~>2)
%@ ;  Matrix = ([]^[0/6]^[2/6]~>1)
%@ ;  Matrix = ([]^[1/6]^[2/6]~>1)
%@ ;  Matrix = ([]^[2/6]^[2/6]~>0)
%@ ;  Matrix = ([]^[3/6]^[2/6]~>0)
%@ ;  Matrix = ([]^[0/6]^[3/6]~>1)
%@ ;  Matrix = ([]^[1/6]^[3/6]~>1)
%@ ;  Matrix = ([]^[2/6]^[3/6]~>0)
%@ ;  Matrix = ([]^[3/6]^[3/6]~>0)
%@ ;  Matrix = ([0/3]^[1/6]^[]~>2)
%@ ;  Matrix = ([]^[0/6]^[2/6]~>1)
%@ ;  Matrix = ([]^[1/6]^[2/6]~>1)
%@ ;  Matrix = ([]^[2/6]^[2/6]~>0)
%@ ;  Matrix = ([]^[3/6]^[2/6]~>0)
%@ ;  Matrix = ([]^[0/6]^[3/6]~>1)
%@ ;  Matrix = ([]^[1/6]^[3/6]~>1)
%@ ;  Matrix = ([]^[2/6]^[3/6]~>0)
%@ ;  Matrix = ([]^[3/6]^[3/6]~>0)
%@ ;  Matrix = ([]^[0/6]^[4/6]~>1)
%@ ;  Matrix = ([]^[1/6]^[4/6]~>1)
%@ ;  Matrix = ([]^[2/6]^[4/6]~>0)
%@ ;  Matrix = ([]^[3/6]^[4/6]~>0)
%@ ;  Matrix = ([]^[0/6]^[2/3]~>1)
%@ ;  Matrix = ([]^[1/6]^[2/3]~>1)
%@ ;  Matrix = ([]^[2/6]^[2/3]~>0)
%@ ;  Matrix = ([]^[3/6]^[2/3]~>0)
%@ ;  Matrix = ([]^[0/6]^[3/3]~>1)
%@ ;  Matrix = ([]^[1/6]^[3/3]~>1)
%@ ;  Matrix = ([]^[2/6]^[3/3]~>0)
%@ ;  Matrix = ([]^[3/6]^[3/3]~>0)
%@ ;  Matrix = ([1/6]^[0/6]^[]~>2)
%@ ;  Matrix = ([1/6]^[1/6]^[]~>2)
%@ ;  Matrix = ([1/6]^[2/6]^[]~>1)
%@ ;  Matrix = ([1/6]^[3/6]^[]~>1)
%@ ;  Matrix = ([1/6]^[1/6]^[]~>2)
%@ ;  Matrix = ([1/6]^[2/6]^[]~>1)
%@ ;  Matrix = ([1/6]^[3/6]^[]~>1)
%@ ;  Matrix = ([1/6]^[4/6]^[]~>1)
%@ ;  Matrix = ([1/6]^[2/3]^[]~>1)
%@ ;  Matrix = ([1/6]^[3/3]^[]~>1)
%@ ;  Matrix = ([]^[2/6,0/0]^[]~>0)
%@ ;  Matrix = ([]^[3/6,0/0]^[]~>0)
%@ ;  Matrix = ([]^[4/6,0/0]^[]~>0)
%@ ;  Matrix = ([]^[2/3,0/0]^[]~>0)
%@ ;  Matrix = ([]^[3/3,0/0]^[]~>0)
%@ ;  false.

%?- J+\(traditional33(T33), D=1, time(findall(M, ccd_d_matrix(T33, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 0.468 seconds
%@    % CPU time: 0.472 seconds
%@    J = 10.

%?- J+\(traditional33(T33), D=2, time(findall(M, ccd_d_matrix(T33, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 2.779 seconds
%@    % CPU time: 2.783 seconds
%@    J = 46.

%?- J+\(traditional33(T33), D=3, time(findall(M, ccd_d_matrix(T33, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 10.030 seconds
%@    % CPU time: 10.034 seconds
%@    J = 154.

%?- J+\(traditional33(T33), D=4, time(findall(M, ccd_d_matrix(T33, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 32.094 seconds
%@    % CPU time: 32.099 seconds
%@    J = 442.

%?- J+\(traditional33(T33), D=5, time(findall(M, ccd_d_matrix(T33, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 90.907 seconds
%@    % CPU time: 90.911 seconds
%@    J = 1162.

%?- J+\(traditional33(T33), D=6, time(findall(M, ccd_d_matrix(T33, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 243.609 seconds
%@    % CPU time: 243.613 seconds
%@    J = 2890.

%?- J+\(traditional33(T33), D=7, time(findall(M, ccd_d_matrix(T33, D, M), Ms)), length(Ms, J)).
%@    % CPU time: 650.724 seconds
%@    % CPU time: 650.728 seconds
%@    J = 6922.
