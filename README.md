# UHF-DIIS
Restricted and Unrestriced Hartree Fock Driver Including the DIIS Algorithm written in Fortran.
The one- and two body parts (h) and (G) tensors were provided by a separate program (TALSH library). 

Test calculation of UHF H2 molecule provided the following output as textfile:



                      ======================================
                               HARTREE-FOCK PROGRAM
                      ======================================

              Input Reader
              -------------------------------------------
              Calculation type:              |  UHF
              Number of electrons:           |  2
              Maximum SCF cycles:            |  15
              SCF convergence threshold:     |  1.000E-13
              Convergence accelerator:       |  DIIS
              Size of the DIIS space:        |  10
              -------------------------------------------



 ----------------------------------------------------------------------------
    Iter         Energy           [F,P]         Error          Conv. Acc.   
 ----------------------------------------------------------------------------
      1        -1.00090511      8.164E-01      1.335E-01
      2        -1.02142176      1.382E-02      2.451E-03          DIIS
      3        -1.02142804      4.252E-06      7.550E-07          DIIS
      4        -1.02142804      1.309E-09      2.324E-10          DIIS
      5        -1.02142804      1.250E-14      5.923E-16          DIIS
 ----------------------------------------------------------------------------

                      *******************************
                                  SUCCES
                       SCF CONVERGED AFTER  5 CYCLES
                      *******************************

           Total SCF Energy                    [a.u.]
           -------------------------------------------------
           Total Energy:              |     -1.0214280393132
           -------------------------------------------------

           Energy Components                   [a.u.]
           -------------------------------------------------
           One electron energy:       |     -1.8910948271135
           Two electron energy:       |      0.4916830675860
           Electronic energy:         |     -1.3994117595275
           Nuclear repulsion:         |      0.3779837202143
           -------------------------------------------------
                            ----------------
                            ORBITAL ENERGIES
                            ----------------

                        NO   OCC          [a.u.]
                        1   1.0000       -0.45386
                        2   1.0000       -0.45386
                        3   0.0000        0.05634
                        4   0.0000        0.05634
                        5   0.0000        0.49611
                        6   0.0000        0.49611
                        7   0.0000        0.58998
                        8   0.0000        0.58998
                        9   0.0000        1.31597
                       10   0.0000        1.31597
                       11   0.0000        1.35871
                       12   0.0000        1.35871
                       13   0.0000        1.35875
                       14   0.0000        1.35875
                       15   0.0000        1.65550
                       16   0.0000        1.65550
                       17   0.0000        1.65554
                       18   0.0000        1.65554

                      *** Terminating Normally ***

