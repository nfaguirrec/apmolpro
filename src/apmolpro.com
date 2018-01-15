!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Copyright (C) 2011-2012 by:                                           !!
!!                                                                          !!
!!    Departamento de Fisica Atomica Molecular y de Agregados               !!
!!    Instituto de Fisica Fundamental                                       !!
!!    Consejo Superior de Investigaciones Cientificas (CSIC)                !!
!!                                                                          !!
!!    Grupo de Quimica Teorica                                              !!
!!    Departamento de Quimica                                               !!
!!    Universidad Nacional de Colombia                                      !!
!!                                                                          !!
!!          Authors:                                                        !!
!!                                                                          !!
!!         - Nestor F. Aguirre                                              !!
!!           nfaguirre@gmail.com                                            !!
!!           nfaguirre@iff.csic.es                                          !!
!!         - Edwin F. Posada                                                !!
!!           efposadac@unal.edu.co                                          !!
!!         - Andres Reyes                                                   !!
!!           areyesv@unal.edu.co                                            !!
!!         - Alexander O. Mitrushchenkov                                    !!
!!           Alexander.Mitrushchenkov@univ-paris-est.fr                     !!
!!         - Maria P. de Lara-Castells                                      !!
!!           pilar.delara.castells@csic.es                                  !!
!!                                                                          !!
!!                            -------------                                 !!
!!                                                                          !!
!!    This program is free software; you can redistribute it and#or modify  !!
!!    it under the terms of the GNU General Public License as published by  !!
!!    the Free Software Foundation; either version 2 of the License, or     !!
!!    (at your option) any later version.                                   !!
!!                                                                          !!
!!    This program is distributed in the hope that it will be useful,       !!
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of        !!
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !!
!!    GNU General Public License for more details.                          !!
!!                                                                          !!
!!    You should have received a copy of the GNU General Public License     !!
!!    along with this program; if not, write to the                         !!
!!    Free Software Foundation, Inc.,                                       !!
!!    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             !!
!!                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

apmolpro={
	show $APMOLPRO_maxit
	show $APMOLPRO_tol
	show $APMOLPRO_hforb
	show $APMOLPRO_dm

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Necesario para inicializar todo
	! lo relativo a HF internamente
	
	{hf
		orbital $APMOLPRO_hforb
	}
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Ejecuta APMO internamente

	APMOLPRO_begin
	
	do iter=1,APMOLPRO_maxit
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Reemplaza H0
		{hf
			orbital $APMOLPRO_hforb
		}

                import 21500.2 JCOUPMTX
                import 21510.2 ICOUPMTX
                import 21520.2 KINMTX

		{matrop
			load Jcoup, SQUARE 21500.2
			load Icoup, SQUARE 21510.2
			load K, SQUARE 21520.2
			add H01, K Jcoup Icoup
			save H01, 21511.2 H0
		}
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Reemplaza la energía nuclear
		
		APMOLPRO_enuc
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se ejecuta un HF para relajar los
		! electrones del sistema con el
		! potencial nuclear incluido
		
		{hf
			save $APMOLPRO_hforb
		}
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! La parte nuclear será mejorada con
		! los electrones a nivel MULTI
		
		APMOLPRO_cMethod
		
		energyMULTI=energy
		
		if( iter.ge.2 ) then
			change = abs(energyMULTI-lastEnergy)
			
			if( change.lt.APMOLPRO_tol ) then
					
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				! En el último cálculo la parte
				! electrónica se lleva hasta
				! nivel CCSD(T)
				
				APMOLPRO_eMethod
				
				text FINAL APMOLPRO ENERGY = $energy
				
				goto endCalc
			endif
			
		endif
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		! Se salvan a archivo los coeficientes
		! (*.vec), la matriz densidad reducida
		! de primer orden (*.dens) y los
		! valores propios de la matriz de fock
		! (*.val) para su posterior uso en APMO
		
		{matrop
			load dens DEN,$APMOLPRO_dm,CHARGE
			load coeff ORB,$APMOLPRO_hforb
			FOCK fock,dens
			tran fockT fock,coeff
			diag evec,eval fockT,1
			save dens 21410.2,SQUARE
			save coeff 21420.2,SQUARE
			save eval 21430.2,VECTOR
		}
		
                export 21410.2 EDEN    ! matrix densidad reducida de primer orden
                export 21420.2 EVEC    ! matrix de coeficientes
                export 21430.2 EVAL    ! energy eigenvalues
		
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Ejecuta APMO internamente
                ! congelando la parte electrónica
                ! para permitir la relajación
                ! nuclear

		APMOLPRO_nrelax
		
		lastEnergy = energyMULTI
	enddo
	
	endCalc
}

