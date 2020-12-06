:- use_module(library(format)).
:- use_module(library(lists)).
:- use_module(library(dcgs)).
:- use_module(library(clpz)).
:- use_module(library(time)).

nth1(N, List, Elt) :- nth0(N, [elt0|List], Elt).
