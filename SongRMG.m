%Coded by Song, S.G. (April 2018)
clear all, close all

disp('###   ')
disp('### Generating Finite Source Models          ###')
disp('###   for Ground Motion Simulation           ###')
disp('### based on 1-Point and 2-Point Statistics  ###')
disp('###   of Kinematic Source Parameters         ###')
disp('###   ')
disp('###        SongRMG, Ver 2018.04              ###')
disp('###   ')
disp('### Author: Seok Goo Song, KIGAM             ###')
disp('### Email : sgsong@kigam.re.kr               ###')


rup.name = 'test1';

% input: fault geometry
rup = gen_src(rup);

% input: target statistics
rup = gen_stats(rup);

%%% End of input parameter generation ~~~~~~~~~~~~~~~~~~

% generate rupture models 
rup = gen_rup(rup);

% generate srf files
%gen_srf(rup);


