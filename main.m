clear all; clc;

% Compounds
Methanol = Compound("Methanol");
Formaldehyde = Compound("Formaldehyde");
Hydrogen = Compound("Hydrogen");
%water = Compound("Water");
%nitrogen = Compound("Nitrogen");

Compounds = [Methanol, Formaldehyde, Hydrogen];

% Reactions
rxn_1 = Reaction(Compounds, [-1, 1, 1]);
rxn_1.CalculateK(873)

Reactions = [rxn_1];

% Reactor
N0 = [1, 1, 0];
P = 1;
T0 = 200 + 273.15; 
Tf = 300 + 273.15;

system = Reactor(Reactions, N0, P, T0, Tf);

























