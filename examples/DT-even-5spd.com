*** H-H
gprint basis
include apmolpro.com

memory 50 M
gthresh twoint=1e-10 prefac=1e-10 energy=1e-10

basis={
	set ORBITAL
	H=aug-cc-pVQZ
	
	set NUCBASIS
s,H,even,nprim=5,ratio=2.5,centre=33.7,dratio=0.8
p,H,even,nprim=5,ratio=2.5,centre=33.7,dratio=0.8
d,H,even,nprim=5,ratio=2.5,centre=33.7,dratio=0.8


	default ORBITAL
}
cartesian

APMOLPRO_maxit = 30
APMOLPRO_tol = 1e-7
APMOLPRO_hforb = 2100.2
APMOLPRO_dm = 21400.2 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First call to APMO in oder
! to build the nuclear wave
! function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_begin={
	{apmo
		species Ha_2,Hb_3
		nucbasis nucbasis,nucbasis
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
		update enuc Ha_2,Hb_3
	}
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Lets to relax the nuclei
! keeping frozen the electrons
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_nrelax={
	{apmo
		load den EDEN
		load orbital EVEC
		load eval EVAL
		frozen e-
		species Ha_2,Hb_3
		nucbasis nucbasis,nucbasis
	}
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Electronic method to use
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_eMethod={
	{apmo
	        molden e-,Ha_2,Hb_3
	}

        {fci
                dm $APMOLPRO_dm
        }
}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nuclear-electron interaction
! method through the
! first-order reduced density
! matrix (record=21400.2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
APMOLPRO_cMethod={
        {fci
                dm $APMOLPRO_dm
        }
}

! Equilibrio BO
!r = 0.74193469

! Optimizacion parcial APMOLPRO
!r = 0.77079904
!r = 0.77081280
r = 0.77091430

symmetry nosym
angstrom
geometry={
        H
        H 1  r
}

{optg procedure=apmolpro gradient=1e-5
}

{property
        density $APMOLPRO_dm
        dm
        qm
}

{matrop
        load D DEN,$APMOLPRO_dm
        natorb nOrb D,1.d-4
        save nOrb 22430.2
}

{put molden DTorb.molden
       orb 22430.2
}

