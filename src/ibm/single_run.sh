#!/usr/bin/env bash 

./xstress 0.01 0.01 0.01 0.01   0.1 0.5 0.01 0.25       0.1 0.1     1.0 0.0 0.0 0.8 0.08 0.5 2.0 1.0 100 100 0.5 1  &

#    mu_feedback  = atof(argv[1]);
#    mu_stress_influx  = atof(argv[2]);
#    mu_influx  = atof(argv[3]);
#    sdmu = atof(argv[4]);
#    s_P_2_NP[0]  = atof(argv[5]);
#    s_P_2_NP[1]  = atof(argv[6]);
#    s_NP_2_P[0]  = atof(argv[7]);
#    s_NP_2_P[1]  = atof(argv[8]);
#    s_12[0]  = atof(argv[9]);
#    s_12[1]  = atof(argv[10]);
#    init_feedback  = atof(argv[11]);
#    init_stress_influx  = atof(argv[12]);
#    init_influx  = atof(argv[13]);
#    cue_P  = atof(argv[14]);
#    cue_NP  = atof(argv[15]);
#    s0  = atof(argv[16]);
#    ad  = atof(argv[17]);
#    aP  = atof(argv[18]);
#    dmax  = atof(argv[19]);
#    zmax  = atof(argv[20]);
#    r = atof(argv[21]);
#    u = atof(argv[22]);
