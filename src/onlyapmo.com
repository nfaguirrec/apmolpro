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
!!           delara@iff.csic.es                                             !!
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

onlyapmo={
	show $APMO_hforb

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Necesario para inicializar todo
	! lo relativo a HF internamente
	
	{hf
		orbital $APMO_hforb
	}
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Ejecuta APMO internamente

	APMO_begin
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Reemplaza H0
	{hf
		orbital $APMO_hforb
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
	
	APMO_enuc
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	! Se ejecuta un HF para relajar los
	! electrones del sistema con el
	! potencial nuclear incluido
	
	{hf
		save $APMO_hforb
	}
}

