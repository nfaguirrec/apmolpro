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

      !**
      ! This is the main subroutine which is directly called
      ! by molpro, when the "apmo" command is found into
      ! the input file
      !
      !     {apmo
      !     }
      !
      !**
      subroutine apmo
        implicit double precision (a-h,o-z)
        include "common/tapes"
        include "common/cgeom"
        include "common/cinput"
        include "common/clib"
        include "common/big"
        include "common/cmpp"
        include "common/thresh"

        character(50) :: baseFileName
        character(5) :: frozenSpecies(ncen+1)
        character(5) :: species(ncen+1)
        character(50) :: nuclearBasis(ncen)
        real(8) :: thresholds(7)
        character(50) :: eBasis
        integer :: ncol, i, j
        character(20) :: token
        logical :: goOut
        logical :: isFrozen
        character(7) :: loadLabel
        character(20) :: loadFileName
        character(5) :: matrixLabel(10)
        character(20) :: matrixOfileName(10)
        integer :: saveMatrixCounter
        integer :: numberOfSpecies
        integer, external :: getpid ! pid of the main molpro process
        character(20) :: idChar

        write(idChar,'(I10)') getpid()
        baseFileName = "input"//trim(adjustl(idChar))

        do j=1,ncen+1
          frozenSpecies(j) = "*NON#"
        end do

        write(iout,'(A)',advance='no')
     >   "1INTERFACE * APMO (Any Particle Molecular Orbital)"
        write(iout,'(7X,A)')
     >   "Authors: N.F. Aguirre, E.F. Posada,"
        write(iout,'(66X,A)',advance='no')
     >   "A. Reyes, M.P. de Lara-Castells"
        write(iout,*) ""

        goOut=.false.
        saveMatrixCounter=1
        numberOfSpecies=0
        species=""
        charge = 0

        thresholds(1)=athr(3) ! 1.0d-12 integralThreshold   TWOINT
        thresholds(2)=athr(20) ! 1.0d-9 electronicDensityMatrixTolerance EDENS
        thresholds(3)=athr(20) ! 1.0d-5 nonElectronicDensityMatrixTolerance EDENS
        thresholds(4)=athr(21) ! 1.0d-9 totalEnergyTolerance  ENERGY
        thresholds(5)=athr(21) ! 1.0d-10 scfElectronicEnergyTolerance  ENERGY
        thresholds(6)=athr(21) ! 1.0d-10 scfNonelectronicEnergyTolerance  ENERGY
        thresholds(7)=athr(1) ! 1.0d-12 doubleZeroThreshold  ZERO

        if(thresholds(2).le.0) thresholds(2)=1.d-6
        if(thresholds(3).le.0) thresholds(3)=1.d-6

        do
          call input(ncol)

          i=1
          do while( i<=ncol )
            call gets(i,token,1)

            select case( token )

              case( 'UPDATE' )
                call gets(i+1,token,1)

                select case( token )
                  case( "ENUC" )
                    write(iout,*) "TYPE PROCESS = UPDATE(ENUC)"
                    do j=1,ncol-1
                      call gets(j+1+i,token,1)
                      numberOfSpecies = numberOfSpecies + 1
                      species(numberOfSpecies) = token
                    end do
                    i=i+ncol

                    if(iprocs .eq. 0) then
                      call updateNuclearEnergy( trim(baseFileName),
     >                 species, numberOfSpecies )
                    end if

                    return
                end select

                i=i+1
                goOut = .true.
                return

              case( 'FROZEN' )
                do j=1,ncol-1
                  call gets(i+j,token,1)
                  frozenSpecies(j) = token
                end do
                i=i+ncol

              case( 'SPECIES' )
                if( ncol-1 /= ncen ) then
                  write (6,'(A,A,I2,A,I2/)', advance='no')
     >              " Number of items in SPECIES",
     >              " should be equal than the number of atoms. ",
     >              ncol-1, "!=", ncen
                  call Error(" ","apmo.f")
                end if

                do j=1,ncen
                  call gets(i+j,token,1)
                  species(j) = token
                end do

                species(ncen+1) = "E-"

                i=i+ncen

              case( 'NUCBASIS' )
                if( ncol-1 /= ncen ) then
                  write (6,'(A,A,I2,A,I2/)', advance='no')
     >              " Number of items in NUCBASIS",
     >              " should be equal than the number of atoms. ",
     >              ncol, "!=", ncen
                  call Error(" ","apmo.f")
                end if

                do j=1,ncen
                  call gets(i+j,token,1)
                  nuclearBasis(j) = token
                end do
                i=i+ncen

              case( 'LOAD' )
                if( ncol-1 /= 2 ) then
                  write (6,'(A,A,I2,A,I2/)', advance='no')
     >              " Number of items in LOAD",
     >              " should be two: label and file",
     >              ncol, "!=", 2
                  call Error(" ","apmo.f")
                end if

                call gets(i+1,token,1)
                loadLabel = token
                call gets(i+2,token,1)
                loadFileName = token

                if(iprocs .eq. 0) then
                  call loadFile( baseFileName, loadLabel, loadFileName )
                end if
                i=i+ncol

              case( 'SAVE' )
                if( ncol-1 /= 2 ) then
                  write (6,'(A,A,I2,A,I2/)', advance='no')
     >              " Number of items in SAVE",
     >              " should be two: matrix and file",
     >              ncol, "!=", 2
                  call Error(" ","apmo.f")
                end if

                call gets(i+1,token,1)
                matrixLabel(saveMatrixCounter) = token
                call gets(i+2,token,1)
                matrixOfileName(saveMatrixCounter) = token

                saveMatrixCounter=saveMatrixCounter+1
                i=i+ncol

              case( 'MOLDEN' )
                write(iout,*) "TYPE PROCESS = MOLDEN"
                do j=1,ncol-1
                  call gets(j+i,token,1)
                  numberOfSpecies = numberOfSpecies + 1
                  species(numberOfSpecies) = token
                  write(iout,*) "Molecular orbitals for "//trim(token)
     >             //" saved"
                end do
                i=i+ncol

                if(iprocs .eq. 0) then
                  call saveMolden( trim(baseFileName),
     >             species, numberOfSpecies )
                end if

                goOut = .true.
                return

              case( 'ENDAPMO' )
                goOut = .true.
                exit

            end select

            i = i+1
          end do

          if( goOut ) exit
        end do

        write(iout,'(A)') ""
        write(iout,'(A)') " THRESHOLDS:"
        write(iout,'(A)') ""
        write(iout,'(A30,ES15.2)')
     >     "INTEGRALS:", thresholds(1)
        write(iout,'(A30,ES15.2)')
     >    "ELECTRONIC DENSITY:", thresholds(2)
        write(iout,'(A30,ES15.2)')
     >    "NON ELECTRONIC DENSITY:", thresholds(3)
        write(iout,'(A30,ES15.2)')
     >    "TOTAL ENERGY:", thresholds(4)
        write(iout,'(A30,ES15.2)')
     >    "SCF ELECTRONIC ENERGY:", thresholds(5)
        write(iout,'(A30,ES15.2)')
     >    "SCF NON ELECTRONIC ENERGY:", thresholds(6)
        write(iout,'(A30,ES15.2)')
     >    "ZERO:", thresholds(7)

        write(iout,'(A)') ""
        write(iout,'(A)') " PARAMETERS:"
        write(iout,'(A)') ""
        write(iout,'(A30,I15)')
     >    "SCF INTERSPECIES MAXIT:", 20000
        write(iout,'(A30,I15)')
     >    "SCF NONELECTRONIC MAXIT:", 20000
        write(iout,'(A30,I15)')
     >    "ITERATION SCHEME:", 3
        write(iout,'(A45)') ""

        call getvar('CHARGE',chrg,unit,ity,nv,1,1)
        write(iout,'(A10,I5)') "CHARGE:", int(chrg)
        write(iout,'(A45)') ""

        write(iout,'(A10,A20,A10,A10)') "SPECIE", "BASIS SET",
     >    "FROZEN", "DUMMY"
        write(iout,'(A)') ""

        do i=1,ncen+1 !! species
          isFrozen=.false.
          do j=1,ncen+1 !! frozenSpecies
              if( species(i) == frozenSpecies(j) ) then
                isFrozen = .true.
                exit
              end if
          end do

          if( species(i) /= "E-" ) then
            if( isFrozen ) then
              write(iout,'(A10,A20,A10)',advance='no') trim(species(i)),
     >          trim(nuclearBasis(i)), "YES"
            else
              write(iout,'(A10,A20,A10)',advance='no') trim(species(i)),
     >          trim(nuclearBasis(i)), "NO"
            end if

            if( charg(i) /= 0.0 ) then
              write(iout,'(A10)') "NO"
            else
              write(iout,'(A10)') "YES"
            end if
          else
            call getstr('BASIS',eBasis,ia,ie,nv,1,1)
            if( isFrozen ) then
              write(iout,'(A10,A20,A10)') "E-", trim(eBasis), "YES"
            else
              write(iout,'(A10,A20,A10)') "E-", trim(eBasis), "NO"
            end if
          end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This code is not parallelized, so
        ! it will be executed only with the
        ! first processor
        if(iprocs .eq. 0) then
          call runAPMO( baseFileName, species,
     >     nuclearBasis, frozenSpecies, thresholds )

          do i=1,(saveMatrixCounter - 1)
           call saveMatrix( baseFileName, matrixLabel(i),
     >      matrixOfileName(i) )
          end do
        end if

        return
      end subroutine apmo

      !**
      ! Runs APMO
      !
      ! @param baseFileName
      !        This is the name of the input file (*.apmo) without extension
      ! @param species
      !        This a vector read from the molpro input file.
      !        For example:
      !           frozen He_4,H,e-
      !        The corresponding vector will be
      !           species = ["HE_4", "H", "E-"]
      !        Warning: The number of species items should
      !        be equal than the number of atoms in the molpro
      !        input file
      ! @param nuclearBasis
      !        This a vector read from the molpro input file.
      !        For example:
      !           nucbasis dirac,2s2p
      !        The corresponding vector will be
      !           species = ["DIRAC", "2S2P"]
      !        Warning: The number of nuclear basis items should
      !        be equal than the number of atoms in the molpro
      !        input file
      ! @param frozenSpecies
      !        This a vector read from the molpro input file.
      !        For example:
      !           frozen He_4,e-
      !        The corresponding vector will be
      !           species = ["HE_4", "E-"]
      !**
      subroutine runAPMO( baseFileName, species,
     >  nuclearBasis, frozenSpecies, thresholds )
        implicit double precision (a-h,o-z)
        include "common/tapes"
        include "common/cgeom"
        include "common/molen"
        include "common/cdirect"

        character(50) :: baseFileName
        character(5) :: frozenSpecies(ncen)
        character(5) :: auxFrozenSpecies(ncen)
        character(5) :: species(ncen)
        integer :: charge
        logical :: showCharge
        character(50) :: nuclearBasis(ncen)
        real(8) :: thresholds(7)
        character(50) :: eBasisFile
        character(50) :: nBasisFile
        character(255) :: apmoLogFile
        character(255) :: apmoDataPath
        character(255) :: path
        character(255) :: pwd
        character(255) :: command
        character(5) :: aux
        character*8 unit
        integer :: numberOfFrozenSpecies
        integer :: i, j, tto
        logical :: isFrozen
        logical :: diracBasis
        integer, external :: getpid ! pid of the main molpro process
        character(20) :: idChar

        write(idChar,'(I10)') getpid()

        !!Basis set name (in lowercase please!)
        eBasisFile="ebasis"//trim(adjustl(idChar))
        nBasisFile="nbasis"//trim(adjustl(idChar))
        apmoLogFile=trim(logfile)//".apmo"

        call getvar('_ANGSTROM',angstrom,unit,ity,nv,1,1)
        call getenv( "MOLPRO_SCR_11", path )
        call getenv( "PWD", pwd )
        call getenv( "APMO_DATA", apmoDataPath )

        call getvar('CHARGE',chrg,unit,ity,nv,1,1)
        charge = int(chrg)

        !!Build apmo input file
        open( unit=34,file=trim(baseFileName)//".apmo",
     >   status='replace', access='sequential', form='formatted' )

        write (34,"(A/,A/,A)")
     >    "SYSTEM_DESCRIPTION='' ",
     >    "",
     >    "GEOMETRY"

        !! Electrons: symbol, basis set and geom
        showCharge=.true.
        do i=1,ncen
          if( charg(i) /= 0.0 ) then
            write(34,"(A7,A20,3F20.15)", advance='no')
     >        "e-("//trim(atname(i))//")", trim(eBasisFile),
     >        rr(:,i)/angstrom

            ! Show electron charge only with the first item
            if( showCharge ) then
              write(34,"(A,I3)") " addParticles=", -charge
              showCharge = .false.
            else
              write(34,'(A)') ""
            end if
          else
            write(34,"(A7,A20,3F20.15)") "e-("//trim(atname(i))//")*",
     >        trim(eBasisFile), rr(:,i)/angstrom
          end if
        end do

        !! Nucleous: symbol, basis set and geom
        do i=1,ncen
          if( nuclearBasis(i)(1:8) == "NUCBASIS" ) then
            if( charg(i) /= 0.0 ) then
              write(34,"(A7,A20,3F20.15)") trim(species(i)),
     >          trim(nBasisFile), rr(:,i)/angstrom
            else
              write(34,"(A7,A20,3F20.15)") trim(species(i))//"*",
     >          trim(nBasisFile), rr(:,i)/angstrom
            end if
          else
            if( charg(i) /= 0.0 ) then
              write(34,"(A7,A20,3F20.15)") trim(species(i)),
     >          trim(nuclearBasis(i)), rr(:,i)/angstrom
            else
              write(34,"(A7,A20,3F20.15)") trim(species(i))//"*",
     >          trim(nuclearBasis(i)), rr(:,i)/angstrom
            end if
          end if
        end do

        write (34,"(A)") "END GEOMETRY"

        write (34,"(A/,A/,A/,A)")
     >    "",
     >    "TASKS",
     >    "      method=RHF",
     >    "END TASKS"

        ! La opcion apmolpro aunque actualmente no esta funcional en APMO
        ! hay que desabilitarla porque presenta problemas
        write (34,"(A/,A/,A,E15.5/,A,E15.5/,A,E15.5/,A,E15.5/,A,E15.5/,
     >    A,E15.5/,A,E15.5/,A,I10/,A,I10/,A/,A/,A/,A)" )
     >    "",
     >    "CONTROL",
     >    "      integralThreshold =",thresholds(1),
     >    "      electronicDensityMatrixTolerance =",thresholds(2),
     >    "      nonElectronicDensityMatrixTolerance =",thresholds(3),
     >    "      totalEnergyTolerance =",thresholds(4),
     >    "      scfElectronicEnergyTolerance=", thresholds(5),
     >    "      scfNonelectronicEnergyTolerance=", thresholds(6),
     >    "      doubleZeroThreshold=", thresholds(7),
     >    "      scfInterspeciesMaximumIterations=", 20000,
     >    "      scfNonelectronicMaxiterations =",20000,
     >    "      transformToCenterOfMass = F",
     >    "      iterationScheme = 3",
     >    "      readCoefficients = T",
     >    "      debugSCFS = T"
!     >    "      apmolpro = T"

        !! Information related with frozen particles
        isFrozen = .FALSE.
        numberOfFrozenSpecies = 0
        auxFrozenSpecies = ""

        do i=1,ncen
          if( frozenSpecies(i) /= "*NON#" ) then
            numberOfFrozenSpecies = numberOfFrozenSpecies + 1
            if(trim(frozenSpecies(i)) == "E-") then
              auxFrozenSpecies(numberOfFrozenSpecies) = "e-"
            else
              aux = ""
              aux = frozenSpecies(i)
              tto = scan( aux, "_")
              do j=2, tto-1
                aux(j:j) = char( IOR( ichar(aux(j:j)),32))
              end do
              auxFrozenSpecies(numberOfFrozenSpecies) = ""
              auxFrozenSpecies(numberOfFrozenSpecies) = aux
            end if
            isFrozen = .TRUE.
          end if
        end do

        if(isFrozen) then
          write(34,"(A,<numberOfFrozenSpecies*2>A)") "      frozen = ",
     >      (trim(auxFrozenSpecies(i))," ", i=1, numberOfFrozenSpecies)
        end if

        write(34,*)
     >    "END CONTROL"

        close(34)

        !! Shows the generated input file
        stat = system("cat "//trim(baseFileName)//".apmo")

        !! Writes the nuclear basis set and
        !! links it to apmo basis library
        nBasisFile = trim(nBasisFile)//".xml"

        diracBasis=.true.
        do i=1,ncen
          if( nuclearBasis(i)(1:8) == "NUCBASIS" ) then
            diracBasis=.false.
            exit
          end if
        end do

        if( .not. diracBasis ) then
          call writeXMLBasis( "NUCBASIS", trim(nBasisFile) )

          stat = system( "repNucBasisSet.sh "//trim(nBasisFile)
     >     //" -ow " )

          stat = system("ln -s $PWD/"//trim(nBasisFile)//" "
     >     //trim(apmoDataPath)//"/basis/ 2> /dev/null")

          stat = system( "cp $PWD/"//trim(nBasisFile)//" "
     >     //trim(pwd)//"/nbasis.xml")
        end if

        !! Writes the electronic basis set and
        !! links it to apmo basis library
        eBasisFile = trim(eBasisFile)//".xml"
        call writeXMLBasis( "ORBITAL", trim(eBasisFile) )

        stat = system("ln -s $PWD/"//trim(eBasisFile)//" "
     >   //trim(apmoDataPath)//"/basis/ 2> /dev/null")

!        stat = system( "cp "//trim(apmoDataPath)//"/basis/"//
!     >   trim(eBasisFile)//" "//trim(pwd)//"/ebasis.xml")

        stat = system( "cp $PWD/"//trim(eBasisFile)//" "
     >   //trim(pwd)//"/ebasis.xml")

        stat = system("apmo -p -i "//trim(baseFileName)//".apmo"//
     >   "  >> "//trim(apmoLogFile))

!        stat = system( "cp $PWD/*.mol "//trim(pwd)//"/" )

        !!Extract and shows the energy components from apmo's output
        command = "awk "//achar(39)//"BEGIN"//achar(123)//"loc=0"//
     >   achar(125)//achar(123)//"if("//achar(36)//"0"//achar(126)//
     >   "/ENERGY COMPONENTS:/)loc=1;if("//achar(36)//"0"//
     >   achar(126)//"/END ENERGY COMPONENTS/)loc=0;if(loc==1)print "//
     >   achar(36)//"0"//achar(125)//"' "//trim(baseFileName)//".out"//
     >   " > outfrag.txt"
        stat = system( command )
        call showFile( "outfrag.txt" )

        !!Extract and shows the dipole moment components from apmo's output
        command = "awk "//achar(39)//"BEGIN"//achar(123)//"loc=0"//
     >   achar(125)//achar(123)//"if("//achar(36)//"0"//achar(126)//
     >   "/DIPOLE:/)loc=1;if("//achar(36)//"0"//
     >   achar(126)//"/END ELECTROSTATIC/)loc=0;if(loc==1)print "//
     >   achar(36)//"0"//achar(125)//"' "//trim(baseFileName)//".out"//
     >   " > outfrag.txt"
        stat = system( command )
        call showFile( "outfrag.txt" )

        !!Extract and shows the expected positions from apmo's output
        command = "awk "//achar(39)//"BEGIN"//achar(123)//"loc=0"//
     >   achar(125)//achar(123)//"if("//achar(36)//"0"//achar(126)//
     >   "/POSITIONS IN ANGSTROMS/)loc=1;if("//achar(36)//"0"//
     >   achar(126)//"/END EXPECTED POSITIONS/)loc=0;if(loc==1)print "//
     >   achar(36)//"0"//achar(125)//"' "//trim(baseFileName)//".out"//
     >   " > outfrag.txt"
        stat = system( command )
        call showFile( "outfrag.txt" )

        stat = system( "rm outfrag.txt" )

        !!remove basis set link in apmo basis set library
        stat = system( "rm "//trim(apmoDataPath)//"/basis/"//
     >   trim(eBasisFile))

        if( nuclearBasis(1)(1:8) == "NUCBASIS" ) then
          stat = system( "rm "//trim(apmoDataPath)//"/basis/"//
     >     trim(nBasisFile))
        end if

      end subroutine runAPMO

      !**
      ! This moves the selected file to the correct name,
      ! selected in the input file of MOLPRO
      !**
      subroutine loadFile( baseFileName, label, ifileName )
        implicit double precision (a-h,o-z)
        include "common/tapes"

        character(*), intent(in) :: baseFileName
        character(*), intent(in) :: label
        character(*), intent(in) :: ifileName

        character(255) :: pwd
        character(255) :: command
        character(255) :: path
        character(20) :: ifileNameEff

        ofileNameEff = ""

        do i=1,len_trim(ifileName)
         ifileNameEff(i:i)=char(IOR(ichar(ifileName(i:i)),32))
        end do

        ifileNameEff( len_trim(ifileName)+1:20 ) = ""

        call getenv( "PWD", pwd )

        select case( trim(label) )

          case( 'ORBITAL' )
            write(iout,*)
     >        " Electronic molecular orbitals readed from "//
     >        trim(pwd)//"/"//trim(ifileNameEff)

            stat = system( "mv "//trim(pwd)//"/"//trim(ifileNameEff)//
     >       " "//trim(baseFileName)//".e-.vec" )

          case( 'DEN' )
            write(iout,*)
     >        " Electronic density readed from "//
     >        trim(pwd)//"/"//trim(ifileNameEff)

            stat = system( "mv "//trim(pwd)//"/"//trim(ifileNameEff)//
     >       " "//trim(baseFileName)//".e-.dens" )

          case( 'EVAL' )
            write(iout,*)
     >        " Electronic molecular orbital eigenvalues readed from "//
     >        trim(pwd)//"/"//trim(ifileNameEff)

            stat = system( "mv "//trim(pwd)//"/"//trim(ifileNameEff)//
     >       " "//trim(baseFileName)//".e-.val" )

          case default
            write(iout,*)
     >        "Unknown option "//trim(label)//" in load"
          call Error(" ","apmo.f")

        end select

        end subroutine loadFile

      !**
      ! This moves the saved matrices by APMO after its execution
      ! to the correct name, selected in the input file of MOLPRO
      !**
      subroutine saveMatrix( baseFileName, matrixLabel, matrixfileName )
        implicit double precision (a-h,o-z)
        include "common/tapes"

        character(*), intent(in) :: baseFileName
        character(*), intent(in) :: matrixLabel
        character(*), intent(in) :: matrixfileName

        character(255) :: pwd
        character(20) :: matrixOfileName

        call getenv( "PWD", pwd )

        matrixOfileName = ""

        do i=1, len_trim(matrixfileName)
         matrixOfileName(i:i)=char(IOR(ichar(matrixfileName(i:i)),32))
        end do

        select case( trim(matrixLabel) )

          case( 'I' )
            write(iout,*)
     >        " Classical-interaction matrix 'I' saved to "
     >        //trim(matrixOfileName)

            stat = system("mv "//trim(baseFileName)//
     >       ".icoup "//trim(matrixOfileName))

          case( 'J' )
            write(iout,*)
     >        " Quantum-interaction matrix 'J' saved to "
     >        //trim(matrixOfileName)

            stat = system("mv "//trim(baseFileName)//".jcoup "//
     >       trim(matrixOfileName))

          case( 'K' )
            write(iout,*)
     >        " Kinetic-energy matrix 'K' saved to "
     >        //trim(matrixOfileName)

            stat = system("mv "//trim(baseFileName)//".kin "//
     >       trim(matrixOfileName))

        end select

        stat = system("cp -rf "//trim(matrixOfileName)//" "//trim(pwd))

        end subroutine saveMatrix

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! This changes to lower case, but keeping the first
        ! letter in upper case. For example: HEA_4 -> Hea_4
        ! it includes E- -> e-
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine toElemName( labelBegin, labelEnd )
          character(5) :: labelBegin
          character(5) :: labelEnd

          character(5) :: buffer
          integer :: i, j, tto

          if( trim(labelBegin) == "E-" ) then
            labelEnd = "e-"
            return
          end if

          buffer = labelBegin
          tto = scan(buffer, "_")

          if( tto == 0 ) then
            tto = len(buffer)+1
          end if

          do j=2, tto-1
            if( buffer(j:j) /= ' ' ) then
!            buffer(j:j) = char( IOR( ichar(buffer(j:j)), 32 ) )
              buffer(j:j) = char( ichar(buffer(j:j))+ 32 )
            end if
          end do

          labelEnd = buffer
        end subroutine toElemName

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Actualiza la energi­a nuclear
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine updateNuclearEnergy( baseFileName, species,
     >     numberOfSpecies )
        implicit double precision (a-h,o-z)
        include "common/tapes"
        include "common/cgeom"
        include "common/molen"
        include "common/cdirect"

        character(*), intent(in) :: baseFileName
        character(5), intent(in) :: species(ncen+1)
        integer, intent(in) :: numberOfSpecies

        character(255) :: path
        character(255) :: pwd
        character(255) :: command
        character(5) :: effSpecies(ncen+1)
        real(8) :: fiener(numberOfSpecies)
        real(8) :: qnkinener(numberOfSpecies)
        real(8) :: Epn(numberOfSpecies)
        real(8) :: coupener(numberOfSpecies)
        character(5) :: buffer
        integer :: i, j, tto
        character*8 unit
        character(5) :: elemName

        write(iout,*) " "

        call getvar('_ANGSTROM',angstrom,unit,ity,nv,1,1)
        call getenv( "MOLPRO_SCR_11", path )
        call getenv( "PWD", pwd )

        do i = 1, numberOfSpecies - 1
          call toElemName( species(i), buffer )
          effSpecies(i) = buffer
        end do

        ! The variables to extract from file
        fiener=0.0_8
        qnkinener=0.0_8
        Epn=0.0_8
        coupener=0.0_8
        ekern=0.0_8

        tcoupener=findString("Coupling energy",
     >   trim(baseFileName)//".out")
        fpotener = findString("Fixed potential energy",
     >   trim(baseFileName)//".out")
        do i = 1, numberOfSpecies - 1
          fiener(i)=findString(trim(effSpecies(i))//
     >     "/Fixed interact. energy=", trim(baseFileName)//".out")

          if( fiener(i) < 1e-10 ) then
            fiener(i)=findString(trim(effSpecies(i))//
     >       "/Fixed interact. energ=", trim(baseFileName)//".out")
          end if

          qnkinener(i)=findString(trim(effSpecies(i))//
     >     " Kinetic energy", trim(baseFileName)//".out")

          Epn(i)=findString(trim(effSpecies(i))//"/"//
     >     trim(effSpecies(i))//" Repulsion energy",
     >     trim(baseFileName)//".out")

          coupener(i)=findString(trim(effSpecies(i))//
     >     " Coupling energy", trim(baseFileName)//".out")
        end do

!        if(numberOfSpecies == 2) then
!          ekern = fiener(1) + fpotener + qnkinener(1) + Epn(1)
!          ekern = fpotener + qnkinener(1) + Epn(1)
!        else
!          ekern = (sum(coupener) - tcoupener) + sum(fiener) + fpotener +
          ekern = sum(coupener) - tcoupener + sum(fiener) + fpotener +
     >     sum(qnkinener) + sum(Epn)
!     >     sum(Epn)
!        end if

        call setvar('ENUC',ekern,'AU',1,1,nv,0)
        call setvar('KNUC',sum(qnkinener),'AU',1,1,nv,0)

        do i = 1, numberOfSpecies -1
          write(iout,"(A40, F20.15)") trim(effSpecies(i))//
     >     " Coupling energy = ", coupener(i)
        end do
        write(iout,"(A40, F20.15)") "- Total coupling energy = ",
     >   -tcoupener
        write(iout,"(A60)") ""

        write(iout,"(A40, F20.15)") "Fixed potential energy = ",
     >   fpotener
        write(iout,"(A60)") ""

        do i = 1, numberOfSpecies -1
          write(iout,"(A40, F20.15)") trim(effSpecies(i))//
     >     "/Fixed interact. energy = ", fiener(i)
        end do
        write(iout,"(A60)") ""

        do i = 1, numberOfSpecies -1
          write(iout,"(A40, F20.15)") trim(effSpecies(i))//
     >     " Kinetic energy = ", qnkinener(i)
        end do
        write(iout,"(A60)") ""

        do i = 1, numberOfSpecies -1
          write(iout,"(A40, F20.15)") trim(effSpecies(i))//"/"//
     >     trim(effSpecies(i))//" Repulsion energy = ", Epn(i)
        end do
        write(iout,"(A60)") ""
        write(iout,"(A60)") ""

        write(iout,"(A40, F20.15)") "Nuclear kinetic"//
     >    " energy to (KNUC) = ", sum(qnkinener)
        write(iout,"(A60)") ""

        write(iout,"(A40, F20.15)") "Nuclear potential"//
     >    " energy to (ENUC) = ", ekern
        write(iout,"(A60)") ""
        write(iout,"(A60)") "*** Warning: ENUC includes KNUC"
        write(iout,"(A60)") ""

          return
        end subroutine updateNuclearEnergy

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Actualiza la energi­a nuclear
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine saveMolden( baseFileName, species,
     >     numberOfSpecies )
        implicit double precision (a-h,o-z)
        include "common/tapes"
        include "common/cgeom"

        character(*), intent(in) :: baseFileName
        character(5), intent(in) :: species(ncen+1)
        integer, intent(in) :: numberOfSpecies

        character(5) :: buffer
        integer :: i

        do i = 1, numberOfSpecies
          call toElemName( species(i), buffer )

          stat = system( "cp $PWD/*."//trim(buffer)//".mol "//
     >     trim(logfile(1:len_trim(logfile)-4))//
     >     "."//trim(buffer)//".molden" )

        end do
      end subroutine saveMolden

      !**
      ! Finds an string in a file and return its value,
      ! following the format "stringToFind = value"
      !**
      function findString(stringToFind, fileName ) result( output )
        implicit none
        character(*), intent(in) :: stringToFind
        character(*), intent(in) :: fileName
        real(8) :: output

        integer :: ssize
        character(255) :: line

        open( UNIT=34,FILE=trim(adjustl(fileName)),
     >   STATUS='unknown',ACCESS='sequential' )

        if (.not.eof(34) ) then
          read(34,'(A)') line
          line=trim( adjustl( trim(line) ) )
        else
          stop "The output file is empty, findString() function"
        end if

        ssize=len_trim( stringToFind )
        do while( line(1:ssize) /= stringToFind(1:ssize)
     >   .and. .not.eof(34) )
          read(34,'(A)') line
          line = trim( adjustl(trim(line)) )
        end do

        if ( line(1:ssize) == stringToFind(1:ssize) ) then
          output=dnum(line(scan(trim(line),"=")+1:len_trim(line)))
        end if

        close(34)
      end function findString

      !**
      ! This subroutine show a file in molpro output file
      !
      ! @param fileName:
      !         The name of the input file to show
      !**
      subroutine showFile( fileName )
        include "common/tapes"
        character(*) :: fileName
        character(255) :: line
        integer :: iostat

        iostat = 1
        open( unit=54, file=fileName, iostat=iostat )
        do while( iostat /= -1 )
          read( 54, "(A)", iostat=iostat ) line
          write( iout, * ) trim(line)
        end do
        close(54)
      end subroutine showFile

      !**
      ! This subroutine extract the species, from
      ! molpro atoms information, for example:
      !
      ! if
      !      atname = ["He", "H", "He", "O" ]
      ! the result will be
      !      species = ["He", "H", "O"]
      !      nSpecies = 3
      !
      ! @output species:
      !         The species without _
      ! @output nSpecies:
      !         The number of species
      !**
      subroutine getAtomTypes( species, nSpecies )
        include "common/cgeom"
        character(5) :: species(ncen)
        integer :: nSpecies

        character(5), allocatable :: temp(:)
        character(5), allocatable :: buffer(:)
        integer :: i, j, k
        logical :: add

        allocate( temp(ncen) )
        allocate( buffer(ncen) )

        species = ""
        buffer = ""
        temp = ""

        k=1
        do i=1,ncen
          add=.true.
          do j=1,ncen
            if(trim(buffer(j)) == trim(atname(i))) then
              add = .false.
            end if
          end do

          if(add) then
            temp(k)=trim(atname(i))
            k = k + 1
          end if

          buffer(i)=trim(atname(i))
        end do

        species(1:k-1) = temp(1:k-1)
        nSpecies = k-1

        deallocate( temp )
        deallocate( buffer )
      end subroutine getAtomTypes

      !**
      ! This subroutine writes the electronic basis
      ! set in xml format
      !
      ! @param basisName:
      !        This is the name basis set
      !
      ! @param fileName:
      !        This is the name of the output file
      !        including the .xml extension
      !**
      subroutine writeXMLBasis( basisName, fileName )
        implicit double precision (a-h,o-z)
        include "common/comun.inc"
        include "common/cgeom"

        character(*) :: basisName
        character(*) :: fileName

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Target variables
        integer :: nAtomTypes !! number of atomTypes
        character(5) :: atomTypes(ncen) !! atom types
        integer :: nContractions(ncen)
        integer :: nPrimitives(ncen,mxfun) !! number of primitives
        character(5) :: shellCode(ncen,mxfun)
        real(8) :: expo(ncen,mxfun)
        real(8) :: coeff(ncen,mxfun)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Auxiliar variables
        character(12) :: terms = "spdfghijklmn"
        logical :: isSpherical
        dimension :: centre(3)
        logical :: procSpecies(ncen)
        integer :: idAtomTypes !! id of atomTypes
        integer :: lastCenter !! last center
        logical :: toShow
        integer :: iContracts !! counter for number of contractions
        integer :: iPrimitives !! counter for primitives
        integer :: i, j, k, l
        real(8) :: expGrp(mxfun)
        real(8) :: coeffGrp(mxfun)
        character(5) :: buffer

        atomTypes = "**"
        procSpecies = .false.
        toShow = .true.

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! It loads the basis set
        write(iout,*) ""
        call load_basis(basisName,basisName,' ',' ',' ',0)
        write(iout,*) ""

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Scanning the basis set from molpro format
        call basis_number( basisName, nbas )
        call getAtomTypes( atomTypes, nAtomTypes )

        call basis_ngroup(basisName,ngrp)
        call basis_spherical(basisName,isSpherical)

        !! It checks that the spherical flag is not active
        if( .not. isSpherical ) then
          write(iout,*) " Using cartesian basis set"
        else
          write(iout,*) " The basis set should be cartesian"
          call Error(" ","apmo.f")
        end if

        lastCenter = 1
        iContracts = 1
        iPrimitives = 1
        idAtomTypes = 1

        do igrp=1,ngrp+1

          if( igrp /= ngrp+1 ) then
            call basis_group(basisName,igrp,nprim,nfun,ioff,
     >        expGrp,centre,icen,minl,maxl,1)
          else
            ! Para que sea diferente a lastCenter y entre
            ! en el siguiente if
            icen = -1
          end if

          !! Evento al momento de alcanzar el cambio de centro
          if( lastCenter /= icen ) then
            do i=1,nAtomTypes
!              if( trim(atomTypes(i)) == trim(atname(lastCenter)) ) then
              if( atname(lastCenter)(1:2) == atomTypes(i)(1:2) ) then

                nContractions(idAtomTypes) = iContracts-1

                iContracts=1
                procSpecies(i)=.true.

                if( toShow ) then
                  idAtomTypes = idAtomTypes+1
                end if

                exit
              end if
            end do

            if( igrp == ngrp+1 ) then
              exit
            end if
          end if

          toShow = .false.
          !! Si la especie no esta procesada se activa para show
          do i=1,nAtomTypes
!            if( trim(atname(icen)) == trim(atomTypes(i)) ) then
            if( atname(icen)(1:2) == atomTypes(i)(1:2) ) then
              if( .not. procSpecies(i) ) then
                toShow = .true.
                exit
              end if
            end if
          end do

          if( toShow ) then
            call basis_contraction_pretty(basisName,igrp,nprim,ncont,
     >       coeffGrp,isSpherical,1)

            do icont=1,ncont
              shellCode(idAtomTypes,iContracts) = terms(minl+1:maxl+1)

              iprimEff=1
              do iprim=1,nprim
                if( abs(coeffGrp((icont-1)*nprim+iprim)) > 1e-6 ) then
                  expo(idAtomTypes,iPrimitives) = expGrp(iprim)

                  ! @todo Esto es temporal, pero por alguna extraña
                  !       razon, las funciona con momento angular
                  !       mayor a 2 (d) las normaliza mal
                  if( nprim == 1 ) then
                    coeff(idAtomTypes,iPrimitives) = 1.0_8
                  else
                    coeff(idAtomTypes,iPrimitives) =
     >               coeffGrp((icont-1)*nprim+iprim)
                  end if

                  iprimEff = iprimEff+1
                  iPrimitives = iPrimitives+1
                end if
              end do

              nPrimitives(idAtomTypes,iContracts) = iprimEff-1
              iContracts = iContracts+1
            end do
          end if

          lastCenter = icen
        enddo

!	write(*,*) "nContractions = ", nContractions
!
!        l=1
!        do i=1,nAtomTypes
!          call toElemName( atomTypes(i), buffer )
!
!          write(*,'(A5,I3)')
!     >     trim(buffer), nContractions(i)
!
!          do j=1,nContractions(i)
!            write(*,'(A,A,I3)') "  ",
!     >        shellCode(i,j), nPrimitives(i,j)
!
!            do k=1,nPrimitives(i,j)
!
!              write(*,'(A,F20.10,F20.10)')
!     >          "    ", expo(i,l), coeff(i,l)
!
!              l=l+1
!            end do
!          end do
!        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Write the basis set into the molpro
        !! output file and the .xml file

        ! @todo Estas lineas fallan para bases nucleares
!         do i=1, len_trim(fileName)
!           fileName(i:i) = char( IOR( ichar( fileName(i:i) ),32))
!         end do

        open( unit=34,file=trim(fileName),
     >    status='replace', access='sequential', form='formatted' )
        write(34,*) "<Basis name="//achar(34)//trim(basisName)//
     >    achar(34)//" type="//achar(34)//"1"//achar(34)//">"

        l=1
        do i=1,nAtomTypes

          call toElemName( atomTypes(i), buffer )

          write(34,'(A)', advance='no') "  <Atom"
          write(34,'(A,A,A)', advance='no') " name="//achar(34),
     >      trim(buffer), achar(34)
          write(34,'(A,A,A)', advance='no') " symbol="//achar(34),
     >      trim(buffer), achar(34)
          write(34,'(A,I3,A)') " numberContractions="//achar(34),
     >      nContractions(i), achar(34)//">"

          write(iout,'(2X,A5,I3)')
     >     trim(buffer), nContractions(i)

          do j=1,nContractions(i)

            write(34,'(A)', advance='no') "    <ContractedGaussian"
            write(34,'(A,I2,A)', advance='no') " size="//achar(34),
     >        nPrimitives(i,j), achar(34)
            write(34,'(A,A,A,I3)') " shellCode="//achar(34),
     >        trim(shellCode(i,j)), achar(34)//">"

            write(iout,'(2X,A,A5,I3)') "  ",
     >        shellCode(i,j), nPrimitives(i,j)

            do k=1,nPrimitives(i,j)

              write(34,'(A)', advance='no') "      <PrimitiveGaussian"
              write(34,'(A,F20.10,A)', advance='no')
     >          " exponent="//achar(34), expo(i,l), achar(34)
              write(34,'(A,F15.10,A)', advance='no')
     >          " coefficient="//achar(34),
     >          coeff(i,l), achar(34)
              write(34,'(A,I4)') "/>"

              write(iout,'(2X,A,F20.10,F15.7)')
     >          "    ", expo(i,l), coeff(i,l)
              l=l+1

            end do

            write(34,'(A)') "    </ContractedGaussian>"
          end do

          write(34,'(A)') "  </Atom>"
        end do

        write(34,*) "</Basis>"
        write(iout,*) ""

        close(34)

        !call system("cat "//trim(fileName))

      end subroutine writeXMLBasis
