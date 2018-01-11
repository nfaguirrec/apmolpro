memory 50 M
include apmolpro.com

APMOLPRO_maxit = 30
APMOLPRO_tol = 1e-6
APMOLPRO_hforb = 2100.2
APMOLPRO_dm = 21400.2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First call to APMO in oder
! to build the nuclear wave
! function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_begin={
	{apmo
		species Ha_3,He_4,He_4
		nucbasis nucbasis,dirac,dirac
		save I,ICOUPMTX
		save J,JCOUPMTX
		save K,KINMTX 
	}
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Updating the nuclear energy
! including the kinetic
! energy from the nuclei
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_enuc={
	{apmo
		update enuc Ha_3
	}
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lets to relax the nuclei
! keeping frozen the electrons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_nrelax={
	{apmo
		load orbital EVEC
		load den EDEN
		load eval EVAL
		frozen e-
		species Ha_3,He_4,He_4
		nucbasis nucbasis,dirac,dirac
	}
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Electronic method to use
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_eMethod={
        {ccsd(t)
        }
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nuclear-electron interaction
! method through the
! first-order reduced density
! matrix (record=21400.2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_cMethod={
        {ccsd
                core 0
                expec relax,dm
                expec dm
                dm $APMOLPRO_dm
                natorb $APMOLPRO_dm
        }
}

basis={
	set ORBITAL

	H=aug-cc-pVQZ
	He=aug-cc-pVQZ

	set NUCBASIS
s,H,even,nprim=5,ratio=2.5,centre=33.7,dratio=0.8
p,H,even,nprim=5,ratio=2.5,centre=33.7,dratio=0.8
d,H,even,nprim=5,ratio=2.5,centre=33.7,dratio=0.8

        s, He, 30.0
        c, 1.1, 1.000000

	default ORBITAL
}
cartesian

r = 0.92491089

set charge=1
symmetry nosym
angstrom
geometry={
         H
        He  1  r
        He  1  r   2  180.0
}

{optg procedure=apmolpro
}

{property
        density $APMOLPRO_dm-10.0
        orbital $APMOLPRO_dm-10.0
        dm
        qm
}

{put molden HeTHeorb.molden
       orb $APMOLPRO_dm-10.0
}

