
 Short instruction how to extract the hadronic (or electromagnetic) 
 shower parameters using CMS the full simulation:

 1) simulate full-sim events:
    -------------------------

    you can chose calorimeter setup for only HCAL or ECAL+HCAL geometry 
    selecting the corresponding line in GenSim_01.py file 

    process.load('Configuration.Geometry.GeometrySimHCAL_cff')
    ##process.load('Configuration.Geometry.GeometrySimECALHCAL_cff')

 2) run job:
    -------- 

    cmsRun GenSim_01.py > test_out01;

 3) extract only hit information from ROOT-file and produce new
    ROOT-file with only shower hit informations: ---------------
    --------------------------------------------

   cmsRun read_01.py > test_out01

 3) there are 3 different examples how to fit the longitudinal shower
    profile: --------------------------------------------------------
    --------

    a) using RooFit package:
       ---------------------
   
       source config_RooFit.csh
       $MyRoot
       .L Pion_RooFit_exam.cxx
        PionFit()

    b) using "standard" fit with Pion_funcFit_exam.C:
       ---------------------------------------------- 

       .L Pion_funcFit_exam.C
        PionFit()

    c) using Minuit package:
       ---------------------
   
       .L fitFunc_loop.C
       .x example_minuit_loop.C

