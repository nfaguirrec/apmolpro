# APMOLPRO

APMOLPRO is an interface between the [APMO](https://sites.google.com/site/lowdinproject/home)
(Any Particle Molecular Orbital) code and the electronic structure package (MOLPRO)[https://www.molpro.net/].
The any particle molecular orbital APMO code
[\[González et al., *Int. J. Quantum Chem.* **108**, 1742 (2008)\]](https://onlinelibrary.wiley.com/doi/abs/10.1002/qua.21584)
implements the model where electrons and light nuclei are treated simultaneously at Hartree-Fock or
second-order Möller-Plesset levels of theory. The APMO -MOLPRO interface allows to include high-
level electronic correlation as implemented in the MOLPRO package and to describe nuclear quantum
effects at Hartree-Fock level of theory with the APMO code.

The examples given In the paper [\[Aguirre et al. *J. Chem. Phys.* **138**, 184113 (2013)\]](http://dx.doi.org/10.1063/1.4803546)
illustrate the use of this implementation on different model systems: <sup>4</sup>
He<sub>2</sub> dimer as a protype of a
weakly bound van der Waals system; isotopomers of \[He–H–He\]<sup>+</sup> molecule as an example of a
hydrogen bonded system; and molecular hydrogen to compare with very accurate non-Born-Oppenheimer calculations.

Flow diagram of the APMO-MOLPRO interface. 
![Terminal](apmolpro.png)

Example of Molpro input file
----------------------------
```
include apmolpro.com
APMOLPRO_maxit = 30
APMOLPRO_tol = 1e-6
APMOLPRO_hforb = 2100.2
APMOLPRO_dm = 21400.2

APMOLPRO_begin={
  {apmo
    species H_1,He_4,He_4
    nucbasis nucbasis,dirac,dirac
    save I,ICOUP
    save J,JCOUP
    save K,KIN
  }
}

APMOLPRO_enuc={
  {apmo
    update enuc H_1
  }
}

APMOLPRO_nrelax={
  {apmo
  load den EDEN
  frozen e-
  species H_1,He_4,He_4
  nucbasis nucbasis,dirac,dirac
  }
}

APMOLPRO_eMethod={
  ccsd(t)
}

APMOLPRO_cMethod={
  {ccsd
    dm $APMOLPRO_dm
  }
}
```

Example of how to specify the electron and nuclear basis set:
-------------------------------------------------------------
```
basis={
  set ORBITAL
  H=aug-cc-pVQZ
  
  set NUCBASIS
  s,H,even,nprim=2,ratio=2.5,centre=33.7,dratio=0.8
  p,H,even,nprim=2,ratio=2.5,centre=33.7,dratio=0.8
  
  default ORBITAL
}
cartesian
```

Example of how to optimize a nuclear basis set:
-----------------------------------------------
```
c1s = 200.0

basis={
  set ORBITAL
  He=aug-cc-pVTZ
  
  set NUCBASIS
  s He c1s
    c 1.1 1.00000
    
  default ORBITAL
}
cartesian

r = 2.99773304

symmetry nosym
angstrom
geometry={
  He
  He 1 r
}

optBasis={
  {optg procedure=apmolpro gradient=1e-5
  }
}

{minimize energy c1s
  method energy simplex,varscale=2,thresh=1e-6,proc=optBasis
}
```
