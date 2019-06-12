%%% Coded by Song, S. (April. 2013)
%%% Set up target 1-point and 2-point statistics

function [rup] = gen_stats(rup)

%%%          [slip Vr Vmax risT];
rup.p1.min = [0   1.0 1.0 0.2];
rup.p1.max = [500 6.0 500 5];

rup.lambda = 1; % min wavelength = 10 km
rup.fN = 6;

rup.p2.tag = 'positive_eigen';

rup.Nstats = 1000;
rup.stats_id = 100;

rup.seed.stats = 1999;

rup.mu = 33e9;

%%% 1-point and 2-point statistics of source parameters

%%%Samping input 1-point and 2-point statistics, given the source statistics model
[rup] = gen_stats_inp(rup);

rup.p1.mu =  [rup.stats(rup.stats_id,1) rup.stats(rup.stats_id,2)*1.5 rup.stats(rup.stats_id,3)];
rup.p1.sig = [rup.stats(rup.stats_id,4) rup.stats(rup.stats_id,5)*1.5 rup.stats(rup.stats_id,6)];


rup.p2.ax = [rup.stats(rup.stats_id,7) rup.stats(rup.stats_id,8)  rup.stats(rup.stats_id,9); ...
             nan                       rup.stats(rup.stats_id,10) rup.stats(rup.stats_id,11); ...
             nan                       nan                        rup.stats(rup.stats_id,12)];

rup.p2.az = [rup.stats(rup.stats_id,13) rup.stats(rup.stats_id,14) rup.stats(rup.stats_id,15); ...
             nan                        rup.stats(rup.stats_id,16) rup.stats(rup.stats_id,17); ...
             nan                        nan                        rup.stats(rup.stats_id,18)];


rup.p2.cc = [1    rup.stats(rup.stats_id,19) rup.stats(rup.stats_id,20); ...
             nan  1                          rup.stats(rup.stats_id,21); ...
             nan  nan                        1];


rup.p2.RDx = [0 0 0; nan 0 0; nan nan 0];
rup.p2.RDz = [0 0 0; nan 0 0; nan nan 0];


