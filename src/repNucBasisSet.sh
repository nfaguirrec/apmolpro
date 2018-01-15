#!/bin/bash
############################################################################
#    Copyright (C) 2011-2012 by:                                           #
#                                                                          #
#    Departamento de Fisica Atomica Molecular y de Agregados               #
#    Instituto de Fisica Fundamental                                       #
#    Consejo Superior de Investigaciones Cientificas (CSIC)                #
#                                                                          #
#    Grupo de Quimica Teorica                                              #
#    Departamento de Quimica                                               #
#    Universidad Nacional de Colombia                                      #
#                                                                          #
#          Authors:                                                        #
#                                                                          #
#         - Nestor F. Aguirre                                              #
#           nfaguirre@gmail.com                                            #
#           nfaguirre@iff.csic.es                                          #
#         - Edwin F. Posada                                                #
#           efposadac@unal.edu.co                                          #
#         - Andres Reyes                                                   #
#           areyesv@unal.edu.co                                            #
#         - Alexander O. Mitrushchenkov                                    #
#           Alexander.Mitrushchenkov@univ-paris-est.fr                     #
#         - Maria P. de Lara-Castells                                      #
#           pilar.delara.castells@csic.es                                  #
#                                                                          #
#                            -------------                                 #
#                                                                          #
#    This program is free software; you can redistribute it and#or modify  #
#    it under the terms of the GNU General Public License as published by  #
#    the Free Software Foundation; either version 2 of the License, or     #
#    (at your option) any later version.                                   #
#                                                                          #
#    This program is distributed in the hope that it will be useful,       #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#    GNU General Public License for more details.                          #
#                                                                          #
#    You should have received a copy of the GNU General Public License     #
#    along with this program; if not, write to the                         #
#    Free Software Foundation, Inc.,                                       #
#    59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
#                                                                          #
############################################################################

IFILE=$1
OVERWRITE=$2
SUFFIXES="a b c d"

getXMLBlock()
{
	local tagXML=$1
	local filter=$2
	local iFile=$3
	
	awk 'BEGIN{loc=0}{ if($0~/<'$tagXML'[[:blank:]]+.*'$filter'.*/) loc=1; if( loc==1 ) print $0; if($0~/<\/'$tagXML'/) loc=0; }' $iFile
}

getXMLLine()
{
	local tagXML=$1
	local filter=$2
	local iFile=$3
	
	awk '($0~/<'$tagXML'[[:blank:]]+.*'$filter'.*>/){ print $0 }' $iFile
}

getIsotopes()
{
	local atomSymbol=$1
	
	local dbFile=$APMO_DATA/dataBases/atomicElements.xml
	
	getXMLBlock "element" "symbol=\"$atomSymbol\"" $dbFile > .isotmp
	grep "isotope" .isotmp | awk 'BEGIN{FS="[[:blank:]=\"]+"}{for(i=1;i<NF;i++) if($i=="massicNumber") printf $(i+1)" " }END{print ""}'
	rm .isotmp
}

atoms=`awk 'BEGIN{FS="[[:blank:]=\"]+"}($0~"<Atom"){ print $6 }' $IFILE`

cat /dev/null > .xmltmp
echo "<!-- Basis generated with repNucBasisSet.sh from APMOLPRO -->" >> .xmltmp
echo "<!-- Author: Nestor Aguirre (nfaguirrec@iff.csic.es) -->" >> .xmltmp
echo "" >> .xmltmp
grep "<Basis " $IFILE >> .xmltmp

for atom in $atoms
do
	isotopes=`getIsotopes $atom $IFILE`

	for A in $isotopes
	do
		getXMLBlock "Atom" "symbol=\"$atom\"" $IFILE | sed '{s/symbol=\"'$atom'\"/symbol=\"'${atom}_$A'\"/g}' >> .xmltmp
		
		for suffix in $SUFFIXES
		do
			getXMLBlock "Atom" "symbol=\"$atom\"" $IFILE | sed '{s/symbol=\"'$atom'\"/symbol=\"'${atom}${suffix}_$A'\"/g}' >> .xmltmp
		done
	done
done

grep "</Basis" $IFILE >> .xmltmp

if [ "$OVERWRITE" = "-ow" ]
then
	mv .xmltmp $IFILE
else
	cat .xmltmp
	rm .xmltmp
fi

