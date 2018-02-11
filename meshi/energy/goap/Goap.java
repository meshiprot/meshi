package meshi.energy.goap;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Scanner;

import meshi.energy.*;
import meshi.molecularElements.Chain;
import meshi.molecularElements.Protein;
import meshi.molecularElements.Residue;
import meshi.molecularElements.ResidueList;
import meshi.molecularElements.atoms.Atom;
import meshi.util.ResidueData;
import meshi.util.Updateable;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;
import meshi.util.info.ChainsInfo;
import meshi.util.info.DoubleInfoElement;
import meshi.util.info.InfoType;

/**
 * The Goap class is a strict (manual, thanks to Tommer Argaman) Java translation of the original Fortran
 * code by Hongyi Zhou and Jeffrey Skolnick Biophysical Journal (2011) 101:2043-2052.
 * All the Fortran code is kept as comments.
 */
public class Goap extends AbstractEnergy implements EvaluatesResidues{
    Protein protein;
    GoapParameters parameters;

    public Goap(Protein protein, GoapParameters parameters, GoapInfo info) {
        super(new Updateable[0],info, EnergyType.NON_DIFFERENTIAL);
        this.protein    = protein;
        this.parameters = parameters;
    }

    public boolean evaluatesResidues() {return true;}


    public EnergyInfoElement evaluate() {
        evaluateResidues(null);
        return info;
    }

    public void evaluateResidues(ChainsInfo chainsInfo) {
        // System.out.println("-------------programStart-----------\n\n\n");



        String filename;
        //char chain;
        double ect_dFIRE, ect2_goapAG, ectTemp, ect2Temp;
        final int maxa_arraySize = 20001, maxres=2001;
        int[] resnum_firstLineOfResInPdbFileIndMinus1 = new int[maxa_arraySize];
        String whole_unknown;
        String resname_residueName; //char*4
        String[] whole1_linesFromFileArray = new String[maxa_arraySize];
        Double prob_errorCharacteristicNumber;
        String ares1_stringResidueNumberFromFile;
        int[][] ibk_idByResPlaceInFileAndAtomPlaceInRes = new int[maxres][16];
        int natom_unknown;
        double[] xp_xCoordinate= new double[maxa_arraySize], yp_yCoordinate= new double[maxa_arraySize], zp_zCoordinate= new double[maxa_arraySize+1];
        Atom[] atoms = new Atom[maxa_arraySize];
        double[] dFireAtom = new double[maxa_arraySize];
        double[] goapAGatom = new double[maxa_arraySize];
        String atmname_atomCode;
        String[] rname_residueNamesByPlaceInFile = new String[maxa_arraySize];
        String[] aname_atomCodesFromFile = new String[maxa_arraySize+1];
        int[] restyp_residueCodeIndicesByPlaceInFile = new int[maxa_arraySize];
        int[] ihflg_flagForNotHydrogen = new int[maxa_arraySize];
        int[] ind1_resCodeIndexByI = new int[maxa_arraySize];
        int[] ind2_atomInResIndexByLine = new int[maxa_arraySize];
        float[] qq_unknown = new float[maxa_arraySize];
        double[][] xn_vectorList=new double[maxa_arraySize][4], xd_vectorList=new double[maxa_arraySize][4];
        double xxh_unknown[]=new double[5],yyh_unknown[]=new double[5],zzh_unknown[]=new double[5];//used before call for dihedral
        int[] ib0_unknown = new int[maxa_arraySize];
        String cct_atomCode2letters;
        double[] xn2_unknown = new double[4], xd2_unknown = new double[4];
        double[] vt_normKtoJ = new double[4], vt2_normJtoK = new double[4];
        int printCounter=1;
//        MeshiWriter finalFeedback;
//        try{finalFeedback = new MeshiWriter("goap_Integrated_Feedback.txt");}
//        catch(Exception ex){throw new RuntimeException(ex);}



//		         do ii=1,nfil
//		c            filename='/home/hzhou2/sparks/1010db/'//afil(ii)
//		            filename=afil(ii) THIS DOESN'T WORK RIGHT. I've changed it in the Fortran version.
//
//		            if(chain.eq.'0') chain=' '
        //chain = ' ';
//		c tommer ect = dFIRE? ect2=goapAG?
//		       ect=0.0
//		       ect2=0.0
        ect_dFIRE=0;
        ect2_goapAG=0;
//		c tommer this opens a pdb structure file
//		 1    open(11,file=filename,status='old') //this is seemingly NOT the start of a loop. 1.4.15
//		          ime=-1
//		c tommer id is set to 1!
//		       id=1
        int id_lineInFileIndex=1;
//		       cx=0
//		       cy=0
//		       cz=0
//		        ires0=0
//		        ires=0
//		        ares0='     '
//		        rnum=0
//		        resnum(1)=0
        double cx_coordinateSum=0,cy_coordinateSum=0,cz_coordinateSum=0;
        int ires0_unknown=0, ires_residueInFileIndex=0, rnum_residuePlaceInFile=0;
        resnum_firstLineOfResInPdbFileIndMinus1[1]=0;
        String ares0_stringResidueInFileNumber = "     ";
//		c tommer it doesn't seem like the variable "line" is used!!! there isn't even such a variable declared.
//		c there is a character variable line1 which is declared but also not used. whatis the meaning of this code line?
//		c does this enable going to the next line in an open file?
//		 15    line = line + 1 //LOOP START!! 1.4.15 "line"????
        for (Atom atom : protein.atoms()) { //debug1
            if (atom.nowhere())continue;
            if (atom.isHydrogen()) continue;
            //while(scanner.hasNextLine()){//LOOP to read from PDB file!
//
//		      read(11,'(a80)',err=10,end=10) whole1(id)
            //err=10 - jumps tp 10 if there is an error on read. 1.4.15
            //whole1_linesFromFileArray[id_lineInFileIndex]=scanner.nextLine();
//		c      if(whole1(id)(22:22).ne.chain.and.ime.eq.1) goto 10 ////////THIS IS THE ONLY PLACE IME_REDUNDANTFLAG IS USED, and it's a comment.
//
//		c tommer to read the next line if doesn't start with 'ATOM'?
//		      if(whole1(id)(1:4).ne.'ATOM') goto 15
//                if (!whole1_linesFromFileArray[id_lineInFileIndex].startsWith("ATOM"))
//                    continue;
//
//
//		c      if(whole1(id)(1:4).ne.'ATOM'.or.whole1(id)(22:22).ne.chain)
//		c     &                     goto 15
//		         read(whole1(id)(18:21),'(a4)') resname
            //resname_residueName=whole1_linesFromFileArray[id_lineInFileIndex].substring(17, 21);
            resname_residueName = atom.residueName();
            //right indices?
//		             iusd=-1
            int iusd_flagForStartingWithResidueCode = -1;

//		          do i=1,20
//		          if(resn(i)(1:3).eq.resname(1:3)) iusd=1
//		          enddo
            for(int residueIndex=1;residueIndex<=20;residueIndex++){
                if(parameters.charges.resn_threeLetterResidueCodes[residueIndex].substring(0, 3).
                        equals(resname_residueName.substring(0, 3)))
                    iusd_flagForStartingWithResidueCode=1;
            }
//		c tommer read the next line if the line is not one of 20 amino acid types
//		          if(iusd.lt.0) goto 15
            if (iusd_flagForStartingWithResidueCode<0) continue;
//
//		      ime=1

//
//		      if(whole1(id)(17:17).ne.' ')
//		     &  READ (whole1(id)(57:60),'(f4.2)') prob   //'(f4.2)' - 4 digits, 2 after decimal
//                if (whole1_linesFromFileArray[id_lineInFileIndex].charAt(16)!=' '){
//                    prob_errorCharacteristicNumber=
//                            Double.parseDouble(whole1_linesFromFileArray[id_lineInFileIndex].substring(56, 60));
            if (!atom.isBackbone()){
                prob_errorCharacteristicNumber=atom.occupancy();

//
//		c tommer checks for reading next line. 1.d-4 - double precision?
//		      if(whole1(id)(17:17).ne.' '.and.prob.lt.0.5)   goto 15
                if (prob_errorCharacteristicNumber<0.5) continue;



//		      if(whole1(id)(17:17).eq.'B'.and.abs(prob-0.5).lt.1.d-4) //1.d-4 -> double precision 1*10^-4 (?)
                //		     &    goto 15
//                    if ((whole1_linesFromFileArray[id_lineInFileIndex].charAt(16)=='B')&&
//                            ((prob_errorCharacteristicNumber-0.5)<0.0001)) //what happens here when we have
//                        // 0.50-0.5?
//                        continue;
            }
//
//		      READ (whole1(id)(23:27),'(a5)') ares1
//                ares1_stringResidueNumberFromFile = whole1_linesFromFileArray[id_lineInFileIndex].substring(22, 27);
            ares1_stringResidueNumberFromFile = String.valueOf(atom.residueNumber());
//		c tommer if read character variable ares1 isn't '     ', fill column ires+1 in array ibk with -1.
//		c tommer ires starts from 0.
//		        if(ares1.ne.ares0) then
//		           ires=ires+1
//		           ares0=ares1
//		           do i=1,15
//		           ibk(ires,i)=-1
//		           enddo
//		        endif
            /*finalFeedback.println("id: "+id_lineInFileIndex+", ares1 :"+ares1_stringResidueNumberFromFile+":, ares0 :"+
            ares0_stringResidueInFileNumber+":\nires: "+ires_residueInFileIndex+ ", atom.residueNumber: "+atom.residueNumber()+
            ", atom.residueName: "+atom.residueName());*/
            if(!ares1_stringResidueNumberFromFile.equals(ares0_stringResidueInFileNumber)){
                ires_residueInFileIndex=ires_residueInFileIndex+1; //ires was initiated with 0.
                ares0_stringResidueInFileNumber=ares1_stringResidueNumberFromFile;
                for (int ibkColIndex=1;ibkColIndex<=15;ibkColIndex++) {
                    ibk_idByResPlaceInFileAndAtomPlaceInRes[ires_residueInFileIndex][ibkColIndex] = -1;
                }
            }// makes sure ares0 and ares1 are equal. 1.4.15
            //distinguishes the switch from one residue to the next (marked by a row of -1's)
//
//
//		c tommer id becomes atom number natom. id starts with 1.
//
//		      natom=id
            natom_unknown=id_lineInFileIndex; //starts with 1
//		c ires becomes residue number rnum.
//		      rnum=ires
            rnum_residuePlaceInFile= ires_residueInFileIndex;
//
//		c tommer value of id is saved in array resnum at (ires+1). why..? from earlier: resnum(1)=0.
//		        resnum(rnum+1)=id
            resnum_firstLineOfResInPdbFileIndMinus1[rnum_residuePlaceInFile+1]=id_lineInFileIndex;
//
//		c           write(*,'(a80)') whole1(id)
//
//		      READ (whole1(id)(31:54),'(3F8.3)') Xp(id),Yp(id),Zp(id)

            xp_xCoordinate[id_lineInFileIndex]= atom.x();
            yp_yCoordinate[id_lineInFileIndex]= atom.y();
            zp_zCoordinate[id_lineInFileIndex]= atom.z();
            atoms[id_lineInFileIndex] = atom;
            //		c tommer coordinates/potential? need to check pdb format.
//		       cx=cx+xp(id)
//		       cy=cy+yp(id)
//		       cz=cz+zp(id)
            cx_coordinateSum=cx_coordinateSum+xp_xCoordinate[id_lineInFileIndex];
            cy_coordinateSum=cy_coordinateSum+yp_yCoordinate[id_lineInFileIndex];
            cz_coordinateSum=cz_coordinateSum+zp_zCoordinate[id_lineInFileIndex];
//
//
//
//
//		         read(whole1(id)(14:17),'(a4)') atmname
//                atmname_atomCode = whole1_linesFromFileArray[id_lineInFileIndex].substring(13, 17);
            atmname_atomCode = atom.name;
            int length = atmname_atomCode.length();
            int diff = 4-length;
            if(diff == 1) { atmname_atomCode += " ";}
            else {
                if (diff == 2) {atmname_atomCode += "  "; }
                else if (diff == 3) atmname_atomCode += "   ";
            }


//		         if(atmname(1:2).eq.'OT') atmname='O   '
            if (atmname_atomCode.substring(0, 2).equals("OT")) atmname_atomCode="O   ";
//		         if(atmname(1:3).eq.'OXT') atmname='O   '
            if (atmname_atomCode.substring(0, 3).equals("OXT")) atmname_atomCode="O   ";
//		         if(resname(1:3).eq.'ILE'.and.atmname(1:3).eq.'CD1')
//		     &           atmname='CD  '
            if (resname_residueName.substring(0, 3).equals("ILE")&&
                    atmname_atomCode.substring(0, 3).equals("CD1"))
                atmname_atomCode="CD  ";

//		         rname(rnum)=resname
            rname_residueNamesByPlaceInFile[rnum_residuePlaceInFile]=resname_residueName;
            //		         aname(id)=atmname
            aname_atomCodesFromFile[id_lineInFileIndex]= atmname_atomCode;
//		          restyp(rnum)=0
            restyp_residueCodeIndicesByPlaceInFile[rnum_residuePlaceInFile]=0; //why??
//
//		c tommer i=id. why is this extra index needed???
//				 i=id

            int i_idDUPLICATEmaybe= id_lineInFileIndex;
//
//
//		c tommer flag for "is not hydrogen"
//		            ihflg(id)=-1
            ihflg_flagForNotHydrogen[id_lineInFileIndex]=-1;
//		            ind2(i)=-1
            ind2_atomInResIndexByLine[i_idDUPLICATEmaybe]=-1;

//		c         if(atmname(1:1).ne.'H'.and.atmname(1:1).ne.'A') then

//		         if(atmname(1:1).ne.'H')then
//		            ihflg(id)=1
            if (atmname_atomCode.charAt(0)!='H'){ //BIG IF! 27.4.15
                ihflg_flagForNotHydrogen[id_lineInFileIndex]=1;

                //		         iok=-1
                int iok_flagIdentifiedAtom=-1;
                //		         do j1=1,20
                //		         if(resn(j1).eq.resname(1:3)) then
                //		c tommer to get from rnum to residue type it is necessary to use restyp(rnum).
                //		            restyp(rnum)=j1
                //		            do j2=1,15
                //		            if(cind(j1,j2).eq.atmname(1:3)) then
                //		c            qq(i)=chg(j1,j2)
                //		            ind1(i)=j1
                //		            ind2(i)=j2
                //		c tommer ibk seems to contain the atom id-s ordered by rnum and place within residue.
                //		            ibk(rnum,j2)=id
                //		c check that non-hydrogen atom is identified ok.
                //		            iok=1
                //		            endif
                //		            enddo
                //		         endif
                //		         enddo
                for(int j1_residueCodeIndex=1;j1_residueCodeIndex<=20;j1_residueCodeIndex++){
                    if (parameters.charges.resn_threeLetterResidueCodes[j1_residueCodeIndex].equals
                            (resname_residueName.substring(0, 3))){
                        restyp_residueCodeIndicesByPlaceInFile[rnum_residuePlaceInFile]=j1_residueCodeIndex;
                        for (int j2_atomInResIndex=1;j2_atomInResIndex<=15;j2_atomInResIndex++){
                            //System.out.println("residueCodeIndex: "+residueCodeIndex+", atomInResInd: "+atomInResIndex);
                            if(parameters.charges.cind_atomsByOrder[j1_residueCodeIndex][j2_atomInResIndex]==null) {
                                throw new RuntimeException("cind entry is null, cannot compare");
                            }
                            //THIS IS TO PREVENT NULL EXCEPTION. this is NOT how it's done in the original code. 27.4.15

                            if (parameters.charges.cind_atomsByOrder[j1_residueCodeIndex][j2_atomInResIndex].equals
                                    (atmname_atomCode.substring(0, 3))){
                                ind1_resCodeIndexByI[i_idDUPLICATEmaybe]=j1_residueCodeIndex;
                                ind2_atomInResIndexByLine[i_idDUPLICATEmaybe]=j2_atomInResIndex;
                                ibk_idByResPlaceInFileAndAtomPlaceInRes[rnum_residuePlaceInFile][j2_atomInResIndex]=id_lineInFileIndex; //here i==id...?

                                //atom in residue is always the same for every residue of the same type. therefore, this array
                                //gives the line in file by the number of the residue in the file and the number of the atom in that residue.
                                iok_flagIdentifiedAtom=1;

                            }
                        }
                    }
                }
                //feedbackWriter.close();

                //
                //		c tommer atom not found. why chose these indices? (10,1)? this looks like a bug. wouldn't disturb the program on normal input...
                //		         if(iok.lt.0) then
                //		c         write(*,*) afil(ii),ires,resname(1:3),atmname(1:3),
                //		c     &            'not found !'
                //		            ind1(i)=10
                //		            ind2(i)=1
                //		            ibk(rnum,1)=id
                //		c         stop
                //		         endif

                if (iok_flagIdentifiedAtom<0){

                    ind1_resCodeIndexByI[i_idDUPLICATEmaybe]=10;// this is 'GLY'?
                    ind2_atomInResIndexByLine[i_idDUPLICATEmaybe]=1;//this is 'N' for Nitrogen??
                    //why GLY Nitrogen?
                    ibk_idByResPlaceInFileAndAtomPlaceInRes[rnum_residuePlaceInFile][1]=id_lineInFileIndex;
                    //why insert the line number in the place of the first atom?
                }
                //
            }
//		         endif
//
//
//		      natom=id
            natom_unknown=id_lineInFileIndex; //isn't this already the case?? [aprox. 120 lines before...]
//		c tommer end of loop to read pdb file
//		      id=id+1
            id_lineInFileIndex++;
//		      goto 15

        } //end of 15 loop
//
//
//		 10   continue
//
//		           close(11)
//
            /*finalFeedback.println("resnum:");
            for (int i=1; i<20;i++)
                finalFeedback.println(i+": "+resnum_firstLineOfResInPdbFileIndMinus1[i]);*/





        int ibb_unknown=0;
//		          do k1=2,rnum
        //ChenDebug loop1
        //new inner loop. goes over residues by order.
        //variables are still read from the pdb file! 4.5.15
        for (int k1_residueInd=2;k1_residueInd<=rnum_residuePlaceInFile;k1_residueInd++){
//		         if((resnum(k1+1)-resnum(k1)).lt.ianum(restyp(k1))) goto 404


            if ((resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd+1]-resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd])<
                    parameters.charges.ianum_numbersOfAtomsInResidues[restyp_residueCodeIndicesByPlaceInFile[k1_residueInd]]) {
//                    break;//checks that there is the right number of getAtoms in the residue. Chen replaced by continue
                continue;
            }


            // ChenDebug loop 2
//		           do k=resnum(k1)+1,resnum(k1+1)
				/*System.out.println("resnum[1]: "+resnum_firstLineOfResInPdbFileIndMinus1[1]+
						"\nresnum[2]: "+resnum_firstLineOfResInPdbFileIndMinus1[2]+
						"\nresnum[3]: "+resnum_firstLineOfResInPdbFileIndMinus1[3]);*/
            for (int k_lineInFileIndex_inResidue=resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd]+1;
                 k_lineInFileIndex_inResidue<=resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd+1];
                 k_lineInFileIndex_inResidue++){
//		           if(ind2(k).gt.0) then



                if(ind2_atomInResIndexByLine[k_lineInFileIndex_inResidue]>0)	//ind2 was given a value
//		           call caldplane(xp,yp,zp,k1,k,ind1(k),ind2(k),ibk,xn1,xd1,ibb)



                    //***************************************now we call CALDPLANE ***********************************************8

//                    if(k_lineInFileIndex_inResidue==10){
//                        finalFeedback.println("\nxp: "+xp_xCoordinate+", "+yp_yCoordinate+", "+zp_zCoordinate+
//                                "\nk1_residueInd: "+k1_residueInd+", k_lineinFilIndex_inResidue: "+k_lineInFileIndex_inResidue+
//                        "\n ind1[k]: "+ind1_resCodeIndexByI[k_lineInFileIndex_inResidue]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInFileIndex_inResidue]+"\n\n");
//                    }

                    ib0_unknown[k_lineInFileIndex_inResidue]=caldplane_unknown_cleaner(xp_xCoordinate,yp_yCoordinate,
                            zp_zCoordinate,
                            k1_residueInd,k_lineInFileIndex_inResidue,ind1_resCodeIndexByI[k_lineInFileIndex_inResidue],
                            ind2_atomInResIndexByLine[k_lineInFileIndex_inResidue],ibk_idByResPlaceInFileAndAtomPlaceInRes);/*,
								Variables.xn1_VectorInput_unknown,Variables.xd1_VectorInput_unknown*///finalFeedback);//ibb is given a value WITHIN THE SUBROUTINE. apparently, unlike in a java function,
                //this value is preserved when leaving the subroutine. Tommer 18.4.15

//		           ib0(k)=ibb  //this is implemented above


                //--------tentative explanations about variables in the function------------
					/*the function gets as input:
					 * 1. coordinate vectors (xp[], yp[], zp[])
					 * 2. place of residue in pdb file [first residue, second residue etc.] (k1)
					 * 3. line-number  of first line in residue (k)
					 * 4. indices of atom and residue code (ind2 and ind1)
					 *
					 *
					 * (internal variables) or RETURNED VARIABLES:
					 *
					 *
					 * (1.) v2: i3-i1, vector from next atom to atom or from previous previous atom to atom
					 * (2. additional notation-not a variable) v4: vector from previous atom atom to atom
					 * 3. xn1: v2+v4 or v4
					 * 4. xd1: (v2Xv4)X(v2+v4) or (xn1Xv2)X xn1
					 *
					 * 5. ibb: the place in the pdb file
					 *
					 *
					 * ADDITIONAL NOTEWORTHY RELATIONS:
					 * if the above descriptions are accurate, it follows that always:
					 * - xn1 and xd1 are orthogonal
					 * - xd1 and (v2+v4) are orthogonal
					 */




                //if (1==1) throw new RuntimeException("this is the end of the checkup");
                //----------------------------------------------------------------------------------------------------------------------
                //------------------------END OF CHECKUP 27.5.15-------------------------------------------//

                //ChenDebug loop3
//			           do l=1,3
                for(int axisInd=1;axisInd<=3;axisInd++){
//			           xn(k,l)=xn1(l)
//			           xd(k,l)=xd1(l)
                    xn_vectorList[k_lineInFileIndex_inResidue][axisInd]=Variables.xn1_Vz_unnormalized[axisInd];
                    xd_vectorList[k_lineInFileIndex_inResidue][axisInd]=Variables.xd1_Vx_unnormalized[axisInd];

//			           enddo   //ChenDebug loop3 ends
                }// create a list of all acquaired xn1, xd1's.
//			           endif
//
            }
//				   enddo  //ChenDebug loop2 ends
        }


//        finalFeedback.println("xn list:");
//            for (int k_ind=1;k_ind<=100;k_ind++) {
//                finalFeedback.println("xn["+k_ind+"]: " + xn_vectorList[k_ind][1]+", "+xn_vectorList[k_ind][2]+", "+xn_vectorList[k_ind][3]);
//            }
//            finalFeedback.println("\n\n");

//		404	continue
//		          enddo
        //ChenDebug loop1 ends
//
        // ChenDebug loop4
//		         do k1=2,rnum-1
        for (int k1_resInd=1; k1_resInd<=rnum_residuePlaceInFile-1;k1_resInd++){      // This line has been changed by Chen 24.8.2017 it was 2 (not 1) I wish I knew why
            //loop over residues, from second in pdb file(i think) to 1 before rnum'th.
//		         if((resnum(k1+1)-resnum(k1)).lt.ianum(restyp(k1))) goto 405
            if((resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd+1]-
                    resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd])<parameters.charges.ianum_numbersOfAtomsInResidues
                    [restyp_residueCodeIndicesByPlaceInFile[k1_resInd]]) {
                // break;  Chen replaced by continue
                continue;
            }
            //check that there's the right number of getAtoms in the residue (bug control?)
            //
            //ChenDebug loop5
//		           do k=resnum(k1)+1,resnum(k1+1)
            for (int k_lineInRes=resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd]+1;k_lineInRes<=
                    resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd+1];k_lineInRes++){
                //loop over the lines of the k1'th residue

//		           cct=aname(k)(1:2)
                cct_atomCode2letters=aname_atomCodesFromFile[k_lineInRes].substring(0, 2);
//		           if(ind2(k).gt.0) then
                if(ind2_atomInResIndexByLine[k_lineInRes]>0){//file not over
                    //if 1.
                    // ChenDebug loop6
//		           do l=1,3
                    for(int l_dimInd=1;l_dimInd<=3;l_dimInd++){
//		           xn1(l)=xn(k,l)
//		           xd1(l)=xd(k,l)
                        Variables.xn1_Vz_unnormalized[l_dimInd]=xn_vectorList[k_lineInRes]
                                [l_dimInd];
                        Variables.xd1_Vx_unnormalized[l_dimInd]=xd_vectorList[k_lineInRes]
                                [l_dimInd];
                    }//regain our local xn1, xd1
//		           enddo
                    //ChenDebug loop6 ends


//
//                   ChenDebug loop7
//		             do j1=k1+1,rnum
                    for (int j1_resCurrentToLast=k1_resInd+1;j1_resCurrentToLast<=//j1key
                            rnum_residuePlaceInFile;j1_resCurrentToLast++)
                    {

//		         if((resnum(j1+1)-resnum(j1)).lt.ianum(restyp(j1))) goto 406
                            /*System.out.println("\n\n\nCheckup start\nj1: "+ j1_resCurrentToLast+", resnum[j1]: "+resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast]+
                            ", ianum[resnum[j1]]: "+parameters.charges.ianum_numbersOfAtomsInResidues[restyp_residueCodeIndicesByPlaceInFile[j1_resCurrentToLast]]+
                            "\nresnum[j1+1]: "+resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast+1]);*/
                        if((resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast+1]-
                                resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast])<parameters.charges.ianum_numbersOfAtomsInResidues
                                [restyp_residueCodeIndicesByPlaceInFile[j1_resCurrentToLast]]) {
                            //break;//zaq4 THIS CAN BREAK THE LOOP AND MAKE TROUBLE. Indeed. Chen replaced by Continue
                            continue;
                        }
//                   ChenDebug loop8
//		             do j=resnum(j1)+1,resnum(j1+1)

                            /*finalFeedback.println("\n\n2nd print, resnum:");
                            for (int i=1; i<20;i++)
                                   finalFeedback.println(i+": "+resnum_firstLineOfResInPdbFileIndMinus1[i]);*/
                            /*if(1==1){
                                finalFeedback.close();
                                throw new RuntimeException("resnum[3]: "+ resnum_firstLineOfResInPdbFileIndMinus1[3]);
                            }*/
                        //if (1==1) throw new RuntimeException("j1: "+j1_resCurrentToLast+
                        //"resnum[j1]: "+resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast]);
                        for(int j_lineInRes=resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast]+1;
                            j_lineInRes<=resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast+1];
                            j_lineInRes++)
                        {
//		             cct=aname(j)(1:2)
//                                cct_atomCode2letters=
//                                        aname_atomCodesFromFile[j_lineInRes].substring(0, 2);
//		           if(ind2(j).gt.0.) then





                                /*if(j_lineInRes==85){
                                    System.out.println("\n\n\nCheckup ind2[j]: " + ind2_atomInResIndexByLine[j_lineInRes]);
                                }*/
                            if(ind2_atomInResIndexByLine[j_lineInRes]>0){//file not over

                                //if 2.

//                   ChenDebug loop9

//		           do l=1,3
                                for(int l_dim=1;l_dim<=3;l_dim++){
//		           xn2(l)=xn(j,l)
//		           xd2(l)=xd(j,l)
                                        /*if(printCounter==1){
                                            finalFeedback.println("\nj="+j_lineInRes+", xn[j-2][l_dim]: "+xn_vectorList[j_lineInRes-2][l_dim]);
                                        }*/
                                    xn2_unknown[l_dim]=xn_vectorList[j_lineInRes][l_dim];
                                    xd2_unknown[l_dim]=xd_vectorList[j_lineInRes][l_dim];
//		           enddo
//                   ChenDebug loop9 ends
                                }


//
//		              rd=sqrt((xp(k)-xp(j))**2+(yp(k)-yp(j))**2+
//		     &                (zp(k)-zp(j))**2)
                                double rd_distanceKJ = Math.sqrt(
                                        (xp_xCoordinate[k_lineInRes]-xp_xCoordinate[j_lineInRes])*(xp_xCoordinate[k_lineInRes]-xp_xCoordinate[j_lineInRes])+
                                                (yp_yCoordinate[k_lineInRes]-yp_yCoordinate[j_lineInRes])*(yp_yCoordinate[k_lineInRes]-yp_yCoordinate[j_lineInRes])+
                                                (zp_zCoordinate[k_lineInRes]-zp_zCoordinate[j_lineInRes])*(zp_zCoordinate[k_lineInRes]-zp_zCoordinate[j_lineInRes]));
                                if (((int)(rd_distanceKJ*2)) >=  parameters.fort21.map_unknown.length) continue;
                                int atomK = k_lineInRes;
                                int atomJ = j_lineInRes;

//
//
//
//		c             if(jj.gt.ibin) jj=ibin
//
//		             jj=map(int(rd*2))
                                int jj_unknown=parameters.fort21.map_unknown[(int)(rd_distanceKJ*2)];//DOES CASTING WORK THE SAME WAY? 27.5
//
//		             if(jj.le.ibin.and.jj.gt.0.1.and.rd.lt.30) then







                                   /* if(j_lineInRes==85){
                                        System.out.println("\nCheckup if 3 - ibin: " + GoapCreator.ibin_unknown+", rd: "+
                                        rd_distanceKJ);
                                    }*/
                                if ((jj_unknown<=GoapCreator.ibin_unknown)&&(jj_unknown>0.1)&&(rd_distanceKJ<30)){

                                    //unknown boubdary check...
                                    //if 3.



//
//		             ee=pot(ind1(k),ind2(k),ind1(j),ind2(j),jj)
//		             ect=ect+ee
                                    double ee_KJjjPotential=parameters.fort21.pot_potential.get_potential(
                                            ind1_resCodeIndexByI[k_lineInRes],
                                            ind2_atomInResIndexByLine[k_lineInRes],
                                            ind1_resCodeIndexByI[j_lineInRes],
                                            ind2_atomInResIndexByLine[j_lineInRes],
                                            jj_unknown);
                                    ectTemp=ect_dFIRE;//checkup 8.6.15
                                    ect_dFIRE=ect_dFIRE+ee_KJjjPotential;
                                    dFireAtom[atomK] += ee_KJjjPotential/2;
                                    dFireAtom[atomJ] += ee_KJjjPotential/2;


//
//		             if((j1-k1).ge.ig) then




                                        /*if(j_lineInRes==85){
                                            System.out.println("\nCheckup if 4 - j1: " + j1_resCurrentToLast+", k1: "+
                                                    k1_resInd+", ig_s: "+parameters.fort31.ig_s_parameter);
                                        }*/
                                    if((j1_resCurrentToLast-k1_resInd)>=parameters.fort31.ig_s_parameter){//only proceed
                                        //if distance between residues is higher or equal to s parameter (7)
                                        //for closer residues angular part is not taken into account
                                        // if 4.
//		             xdd=1./rd
                                        double xdd_distReciprocal= 1./rd_distanceKJ;
//		             vt(1)=(xp(j)-xp(k))*xdd   ! 1->2
//		             vt(2)=(yp(j)-yp(k))*xdd   ! 1->2
//		             vt(3)=(zp(j)-zp(k))*xdd   ! 1->2
                                        vt_normKtoJ[1] = (xp_xCoordinate[j_lineInRes]-
                                                xp_xCoordinate[k_lineInRes])*xdd_distReciprocal;
                                        vt_normKtoJ[2]=(yp_yCoordinate[j_lineInRes]-
                                                yp_yCoordinate[k_lineInRes])*xdd_distReciprocal;
                                        vt_normKtoJ[3]=(zp_zCoordinate[j_lineInRes]-
                                                zp_zCoordinate[k_lineInRes])*xdd_distReciprocal;

//
//		             vt2(1)=-vt(1)
//		             vt2(2)=-vt(2)
//		             vt2(3)=-vt(3)
                                        vt2_normJtoK[1]=-vt_normKtoJ[1];
                                        vt2_normJtoK[2]=-vt_normKtoJ[2];
                                        vt2_normJtoK[3]=-vt_normKtoJ[3];




                                        //feedbackWriter.close();
                                        //if (1==1) throw new RuntimeException("this is the end of the checkup 27.5");
//

                                        //-----------------------------------------------------------27.5.15-------------------------------------------------------------------

//		             call calang2(xn1,vt,ang1,cs1)
                                        double  pointerToAngAndCos1_thetaA_cosThetaA[]=new double[3], pointerToAngAndCos2_thetaB_cosThetaB[]=new double[3];
                                        int mm1_cosThetaA_xn1_vt_bin, mm2_psyA_xd1_vt_bin, mm3_cosThetaB_xn2_vt2_bin, mm4_psyB_xd2_vt2_bin;

                                        calang2_calculateAngleInDeg(Variables.xn1_Vz_unnormalized,vt_normKtoJ, pointerToAngAndCos1_thetaA_cosThetaA);


//		             mm1=int((cs1+1.001)*6.0)+1
                                        mm1_cosThetaA_xn1_vt_bin= (int)((pointerToAngAndCos1_thetaA_cosThetaA[2]+1.001)*6.0)+1;
//		             if(mm1.gt.12) mm1=12
                                        if(mm1_cosThetaA_xn1_vt_bin>12) mm1_cosThetaA_xn1_vt_bin=12;
//		             call calphi(xn1,xd1,vt,phi1)
                                        double phi1_signedAngleBet_xd1Andvt_psyA=
                                                calphi_signedAngleBet_xh1ANDxd_PSY(Variables.xn1_Vz_unnormalized,
                                                        Variables.xd1_Vx_unnormalized,vt_normKtoJ);

//		             mm2=int((phi1+180.001)/30.0)+1
                                        mm2_psyA_xd1_vt_bin= (int)((phi1_signedAngleBet_xd1Andvt_psyA+
                                                180.001)/30.0)+1;
//		             if(mm2.gt.12) mm2=12
                                        if(mm2_psyA_xd1_vt_bin>12) mm2_psyA_xd1_vt_bin=12;
//----------------------------------------------------------------1.6.15 checkup----------------------------------------




                                        //feedbackWriter.close();
                                        //if (1==1) throw new RuntimeException("this is the end of the checkup 1.6");




//		             call calang2(xn2,vt2,ang2,cs2)
//                                            if(printCounter==1){
//                                                finalFeedback.println("\nxn2: "+xn2_unknown[1]+", "+xn2_unknown[2]+", "+xn2_unknown[3]+"\nvt2: "+
//                                                        vt2_normJtoK[1]+", "+vt2_normJtoK[2]+", "+vt2_normJtoK[3]+"\n");
//                                            }
                                        calang2_calculateAngleInDeg(
                                                xn2_unknown, vt2_normJtoK, pointerToAngAndCos2_thetaB_cosThetaB);




//		             mm3=int((cs2+1.001)*6.0)+1
                                        mm3_cosThetaB_xn2_vt2_bin=(int)((pointerToAngAndCos2_thetaB_cosThetaB[2]+1.001)*6.0)+1;
//		             if(mm3.gt.12) mm3=12
                                        if(mm3_cosThetaB_xn2_vt2_bin>12) mm3_cosThetaB_xn2_vt2_bin=12;
//		             call calphi(xn2,xd2,vt2,phi2)
                                        double phi2_signedAngleBet_xd2Andvt2_psyB=calphi_signedAngleBet_xh1ANDxd_PSY(
                                                xn2_unknown, xd2_unknown, vt2_normJtoK);
//		             mm4=int((phi2+180.001)/30.0)+1
                                        mm4_psyB_xd2_vt2_bin=(int)((phi2_signedAngleBet_xd2Andvt2_psyB+180.001)/30.0)+1;

//		             if(mm4.gt.12) mm4=12
                                        if(mm4_psyB_xd2_vt2_bin>12) mm4_psyB_xd2_vt2_bin=12;
//
//--------------------------------------------------1.6 checkup---------------------------------------




//			   		xxh(1)=xp(j)+xn2(1)
//		             yyh(1)=yp(j)+xn2(2)
//		             zzh(1)=zp(j)+xn2(3)
//		             xxh(2)=xp(j)
//		             yyh(2)=yp(j)
//		             zzh(2)=zp(j)
//		             xxh(3)=xp(k)
//		             yyh(3)=yp(k)
//		             zzh(3)=zp(k)
//		             xxh(4)=xp(k)+xn1(1)
//		             yyh(4)=yp(k)+xn1(2)
//		             zzh(4)=zp(k)+xn1(3)
                                        xxh_unknown[1]=xp_xCoordinate[j_lineInRes]+xn2_unknown[1];
                                        yyh_unknown[1]=yp_yCoordinate[j_lineInRes]+xn2_unknown[2];
                                        zzh_unknown[1]=zp_zCoordinate[j_lineInRes]+xn2_unknown[3];
                                        xxh_unknown[2]=xp_xCoordinate[j_lineInRes];
                                        yyh_unknown[2]=yp_yCoordinate[j_lineInRes];
                                        zzh_unknown[2]=zp_zCoordinate[j_lineInRes];
                                        xxh_unknown[3]=xp_xCoordinate[k_lineInRes];
                                        yyh_unknown[3]=yp_yCoordinate[k_lineInRes];
                                        zzh_unknown[3]=zp_zCoordinate[k_lineInRes];
                                        xxh_unknown[4]=xp_xCoordinate[k_lineInRes]+Variables.xn1_Vz_unnormalized[1];
                                        yyh_unknown[4]=yp_yCoordinate[k_lineInRes]+Variables.xn1_Vz_unnormalized[2];
                                        zzh_unknown[4]=zp_zCoordinate[k_lineInRes]+Variables.xn1_Vz_unnormalized[3];
//



//		            call dihedral(xxh,yyh,zzh,dih)
                                        double dih_angleBetweenNormsToConsecutivePlanes_xhi=
                                                dihedral_calc_Xhi(xxh_unknown,yyh_unknown,zzh_unknown);








											/*if(printCounter==11172){//checkup 30.6.15
												System.out.println("checkup 30.6\n"+"dih: "+dih_angleBetweenNormsToConsecutivePlanes);
												System.out.println("xxh: "+xxh_unknown[1]+", "+xxh_unknown[2]+", "+xxh_unknown[3]+", "+xxh_unknown[4]);
												System.out.println("yyh: "+yyh_unknown[1]+", "+yyh_unknown[2]+", "+yyh_unknown[3]+", "+yyh_unknown[4]);
												System.out.println("zzh: "+zzh_unknown[1]+", "+zzh_unknown[2]+", "+zzh_unknown[3]+", "+zzh_unknown[4]);
												System.out.println("\n");
											}*/



//		             mm5=min(int((dih+180.001)/30.0)+1,12)
                                        int mm5_dih_xhi_bin=Math.min(
                                                (int)((dih_angleBetweenNormsToConsecutivePlanes_xhi+180.001)/30.0)+1, 12);
//
//
//
//
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm1,1)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm2,2)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm3,3)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm4,4)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm5,5)
                                        ect2Temp=ect2_goapAG;
                                        double eGoapAG;
                                        if(printCounter==1) {
                                            eGoapAG = parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm1_cosThetaA_xn1_vt_bin, 1);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;
//                                                finalFeedback.println("\nfirst addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp));

                                            ect2Temp=ect2_goapAG;
                                            eGoapAG =  parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm2_psyA_xd1_vt_bin, 2);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;
//                                                finalFeedback.println("\nsecond addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp));

                                            ect2Temp=ect2_goapAG;
                                            eGoapAG = parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm3_cosThetaB_xn2_vt2_bin, 3);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;
//                                               finalFeedback.println("\nthird addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                "\nk: "+k_lineInRes+", j: "+j_lineInRes+
//                                                "\nind1[k]: "+ind1_resCodeIndexByI[k_lineInRes]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInRes]+
//                                                "\nind1[j]: "+ind1_resCodeIndexByI[j_lineInRes]+", ind2[j]: "+ind2_atomInResIndexByLine[j_lineInRes]+
//                                                "\njj: "+jj_unknown+", mm3: "+mm3_cosThetaB_xn2_vt2_bin);

                                            ect2Temp=ect2_goapAG;
                                            eGoapAG =  parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm4_psyB_xd2_vt2_bin, 4);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;
//                                                finalFeedback.println("\nfourth addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                "\nk: "+k_lineInRes+", j: "+j_lineInRes+
//                                                "\nind1[k]: "+ind1_resCodeIndexByI[k_lineInRes]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInRes]+
//                                                "\nind1[j]: "+ind1_resCodeIndexByI[j_lineInRes]+", ind2[j]: "+ind2_atomInResIndexByLine[j_lineInRes]+
//                                                "\njj: "+jj_unknown+", mm4: "+mm4_psyB_xd2_vt2_bin);
                                            ect2Temp=ect2_goapAG;
                                            eGoapAG =  parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm5_dih_xhi_bin, 5);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;
//                                                finalFeedback.println("\nfifth addition"+
//												"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                "\nk: "+k_lineInRes+", j: "+j_lineInRes+
//                                                "\nind1[k]: "+ind1_resCodeIndexByI[k_lineInRes]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInRes]+
//                                                "\nind1[j]: "+ind1_resCodeIndexByI[j_lineInRes]+", ind2[j]: "+ind2_atomInResIndexByLine[j_lineInRes]+
//                                                "\njj: "+jj_unknown+", mm5: "+mm5_dih_xhi_bin);


                                        }
                                        else{
                                            eGoapAG = parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm1_cosThetaA_xn1_vt_bin, 1);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;

                                            eGoapAG = parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm2_psyA_xd1_vt_bin, 2);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;

                                            eGoapAG = parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm3_cosThetaB_xn2_vt2_bin, 3);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;

                                            eGoapAG = parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm4_psyB_xd2_vt2_bin, 4);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;

                                            eGoapAG = parameters.fort31.cnttheta_unknown.get_cnttheta(
                                                    ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                    ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                    jj_unknown, mm5_dih_xhi_bin, 5);
                                            ect2_goapAG += eGoapAG;
                                            goapAGatom[atomK] +=eGoapAG/2;
                                            goapAGatom[atomJ] +=eGoapAG/2;
//                                                finalFeedback.println("\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                    "\ncount: "+printCounter);
                                        }


											/*finalFeedback.println("\nect: "+ect_dFIRE+
													"\nect increment: "+(ect_dFIRE-ectTemp)+
													"\nect2: "+ect2_goapAG+
													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
													"\ncount: "+printCounter);*/

                                        printCounter++;




                                        //feedbackWriter.close();
                                        //if (1==1) throw new RuntimeException("this is the end of the checkup 3.6");


//		             endif ! .ge.ig
                                    }//if (j1-k1)> ig
//
//		             endif
                                }// boundary check for jj and rd
//
//		            endif
                            }
//		1234              continue
//		            enddo
//                   ChenDebug loop8ends



                        }
                        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@zaq3@@@@@@@@@@@
                        //this is the last place where the fortran program
                        //can be exited without error. 6.5.15


//		406	continue
//		            enddo
//                  ChenDebug loop7 ends
                    }
//
//		           endif
                }
//		          enddo
//                  ChenDebug loop5 ends
            }
//
//
//		405	continue
//		         enddo
//                  ChenDebug loop4 ends
        }
//
//		       write(*,555) ii,' ',afil(ii),ect+ect2,ect,ect2
//		403	continue
//		      enddo
//		c tommer this is the end of the main loop??
//
//		555	format(1x,i5,a1,1x,a30,1x,f12.2,1x,2f12.2)
//

        //System.out.println("\n\n\nOutput, 30.6.15:\nii: "+ii_fileNumber+", afil[ii]: "+afil_filenames[ii_fileNumber]+
        //        "\ngoap score: "+(ect_dFIRE+ect2_goapAG)+"\ndFIRE score: "+ect_dFIRE+"\ngoapAG score: "+ect2_goapAG);//checkup
        //finalFeedback.println("\n\n\nOutput, 6.6.15:\nii: "+ii_fileNumber+", afil[ii]: "+afil_filenames[ii_fileNumber]+
        //	"\ngoap score: "+(ect_dFIRE+ect2_goapAG)+"\ndFIRE score: "+ect_dFIRE+"\ngoapAG score: "+ect2_goapAG);//checkup 8.6.15




//        stop
//      end

//        finalFeedback.close();



        info.setValue(-ect_dFIRE-ect2_goapAG);
        ((GoapInfo) info).dFIRE.setValue(-ect_dFIRE);
        ((GoapInfo) info).goapAG.setValue(-ect2_goapAG);
        if (chainsInfo != null) {
            ResidueData dFireData = new ResidueData(chainsInfo);
            ResidueData goapAGdata = new ResidueData(chainsInfo);
            for (int iAtom = 0; iAtom < maxa_arraySize; iAtom++) {
                if (atoms[iAtom] != null) {
                    Atom atom = atoms[iAtom];
                    dFireData.add(atom,dFireAtom[iAtom]);
                    goapAGdata.add(atom,goapAGatom[iAtom]);
                }
            }

            for (int iChain = 0; iChain < protein.chains().size(); iChain++) {
                Chain chain = protein.chains().get(iChain);
                for (int jResidue = 0; jResidue < chain.size(); jResidue++) {
                    Residue residue = chain.get(jResidue);
                    if (!residue.dummy()) {
                        chainsInfo.get(iChain).get(jResidue).add(
                                new DoubleInfoElement(InfoType.D_FIRE, "dFire",   dFireData.get(iChain, jResidue)/residue.getAtoms().size()));
                        chainsInfo.get(iChain).get(jResidue).add(
                                new DoubleInfoElement(InfoType.GOAP_AG, "goapAG",   goapAGdata.get(iChain, jResidue)/residue.getAtoms().size()));
                    }
                }
            }
        }
    }

    public void evaluateAtoms(){

    }
    public void test(TotalEnergy energy, Atom atom) {}



    /**
     * A standalone version.
     * @param args  - command line arguments.
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

//		cfil  The output  : 	ect is DFIRE, ect2 is GOAP_AG
//		c              GOAP=DFIRE+GOAP_AG
//   	c	all inputs are at /gpfs1/active/hzhou2/fold/dplanedfire/ 
//		c       which was set in variable base
//		c	 basic flow:  (1) read in structure list 
//		c                     (2) read in potential data
//		c	              (3) loop over each structure
//		c	       inside the loop: (a) read in structure
//		c	                        (b) calculate local system vectors for each atom
//		c	                        (c) loop over all pairs and evaluate
//		CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
//		      implicit real*8(a-h,o-z)
//		      integer maxa,nwt,newwt,maxp,natom,numdel,
//		     &        ntyp,maxit,maxres
//		      logical iwrit,ioutwt,slmut,slmut1
//		      real*8 cbox,cx,cy,cz,vdwradi(5),pi,pi3,theta,rcut,sgcut,rd
//		      parameter(maxa=20000,cbox=18.856,maxit=100,
//		     &          ntyp=20,natyp=8,maxres=1000,ibin=20)
//		      integer idel(maxa),resnum(maxa),restyp(maxa),pol(20),
//		     &        ihflg(maxa),rnum,ncel,iused(ntyp),imut(maxres),
//		     &        itgt(maxres),iatyp(maxa),iu(natyp,natyp),
//		     &        mutmap(ntyp,12),ind1(maxa),ind2(maxa),
//		     &        map(500),ibk(maxres,15),ianum(20),ib0(maxa)
//		      real*8 xp(maxa),yp(maxa),zp(maxa),ddg(maxres),resdep(maxres),  
//		     &        hd(20),asa(20),qq(maxa),rco(maxres),ele(maxa),aa
//		      real*8 xpn(maxa),ypn(maxa),zpn(maxa),atmdep(maxa),hvnum(ntyp),
//		     &       icnt(50),pot(20,15,20,15,ibin),rrme(ibin),rghist(1000),
//		     &       scalco(ibin),ddr(ibin),pb(20,15,20,15,ibin)
//		      real*8 sx,sy,sz,tx,ty,tz,xm(3,4),dep(maxa,maxit),avedep,
//		     &       entropy(ntyp),volume(ntyp),aenv(maxres),para(ntyp),
//		     &       para1(ntyp),alf(natyp,natyp),bkbn(natyp,natyp),
//		     &       chg(ntyp,15),rmsd(3000),cntnum(20,15,20,15,ibin),
//		     &       ttlcnt(ibin),xmol(20,15),expnum(20,15,20,15,ibin),
//		     &       unk1(ibin,ibin),unk2(ibin),hnum,ha1,ha2,dh,hnumex,
//		     &       ymol(20,15),pbtl(ibin),xn(maxa,3),xd(maxa,3),
//		     &       xn2(3),xd2(3),xn1(3),xd1(3),cnttheta(20,15,20,15,ibin,12,5)
//		     &      ,vt(3),vt2(3),cnttheta0(12,5),av0(5),sd0(5),ppx(12)
//		     &      ,xxh(4),yyh(4),zzh(4),dih,dih2
//		      character resname*4,atmname*4,ctmp*4,resn(ntyp)*3,rname(maxa)*4,
//		     &    aname(maxa)*4,line1*80,cct*2,chain*1,ares0*5,ares1*5
//		     &  ,ctmp1*4,ctmp2*2,ctmp3*4
//		       character resnm*3,atnm*3,cind(ntyp,20)*3,base*70
//
//		       data qq/maxa*0.0/,ele/maxa*0.0/,aa/0.0/
//		       data map/500*0/
//		       data cind/400*'XOP'/
//		       common /fff/ips,base
//
//		       data asa/360.,425.,429.,402.,399.,375.,470.,447.,347.,
//		     &  323,370.,352.,408.,379.,405.,378.,432.,468.,445.,364./
//
//		       data hd/4.19,4.80,6.03,5.87,5.89,4.51,7.14,4.90,2.03,1.01,
//		     &     2.53,1.59,2.48,1.40,1.78,1.09,2.95,1.78,2.01,3.51/
//
//
//
//		      data hvnum/2,4,7,4,4,3,10,8,1,0,
//		     &           3,2,5,4,5,4,6,7,5,3/               
//
//		      data pol/1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,0/
//
//		      data volume/45.,104.,130.,101.,101.,75.,168.,133.,26.,0.,
//		     &     56.,30.0,86.0,64.,77.,53.,96.0,129.,106.,59./
//
//
//
//		      data entropy/-0.55,-1.61,-0.58,-0.89,-0.78,-0.51,-0.97,-0.98,
//		     & 0.0,0.0,-1.63,-1.71,-2.11,-1.57,-1.80,-1.25,-0.96,-2.03,
//		     & -1.94,0.0/            
//
//		      data vdwradi,pi/2.0,1.85,1.70,2.0,0.,3.1415926/
//
//		      data resn/'CYS','MET','PHE','ILE','LEU','VAL','TRP','TYR',
//		     &          'ALA','GLY','THR','SER','GLN','ASN','GLU','ASP',
//		     &          'HIS','ARG','LYS','PRO'/
//		   
//
//		      common /ddd/icnt
//
//
        MeshiLineReader reader_input = new MeshiLineReader(".\\tmp.inp.txt");
//		MeshiWriter finalFeedback = new MeshiWriter("..\\goap_Main_Feedback.txt");
//		      character*80 whole, filename*100,whole1(maxa),fil1*40,fil2*40,
//		     &              afil(60000)*80
//
//		c chen base is the directory of the goap input data, see readme.
//
        //Variables.base_goapalonePath = reader_input.readLine();
        Variables.base_goapalonePath = "C:\\Users\\Tomer\\Program 2015\\meshi9.10\\tests\\goapTest\\t0281";
        //this is not obviously necessary... do we need the path? is this the way to get it? Tommer 3.2.15
//		          read(*,'(a70)') base
//
//		c chen a loop that reads an input file into an array of characters "afil", input file names, see readme.
//		        i=1 
//		c!!!        open(unit=20,file='pdbfile.list',status='old')
//		 2001   read(*,*,end=2000) afil(i)
//		        i=i+1
//		        goto 2001
//		c chen done reading the file names.
        String[] afil_filenames= new String[6000];
        int fileNumberIndex;
        for (fileNumberIndex=1; reader_input.ready(); fileNumberIndex++){
            afil_filenames[fileNumberIndex]= reader_input.readLine();
        }
        reader_input.close();

        int nfil_numberOfFiles = fileNumberIndex-1;




//		 2000   continue
//		c chen nfil is the number of file names
//		        nfil=i-1
//		c        close(20)
//---------------------------------------------------------------------------------2
//
//		            do i=1,500
//		            map(i)=0
//		            enddo
        int[] map_unknown = new int[501]; //this should be initializied with zeros by default. 11.3.15
        //the array should start from map[1] hence the size should be 501.
        //will later be filled (to 50'th slot) by numbers from beginning of "fort" files.

//		            hcut=3.5 // initialization without declaration for int. interesting fortran feature.
        double hcut_unknown = 3.5;

//		            dh=2./ibin
        final int ibin_unknown=20; //this is originaly declared inside a "pararmeter(...)"
        double dh_unknown = 2./ibin_unknown;


        int ips_goapPathLength = Variables.base_goapalonePath.length();




//
//		c         base='/home/hzhou2/fold/dplanedfire '
//		c         base='/gpfs1/active/hzhou2/fold/dplanedfire '
//		        ips=index(base,' ')-1
//		c tommer ips = index of first ' ' in base (-1)=length of the path to the goap-alone directory
//
//		c        open(unit=20,file=base(1:ips)//'/fort.21_1.61_10',
//
//		        open(unit=20,file=base(1:ips)//'/fort.21_1.61_2',

        MeshiLineReader data_from_fort21Reader_20 = new MeshiLineReader(Variables.base_goapalonePath+
                "\\fort.21_1.61_2");
        data_from_fort21Reader_20.readLine();//read comment line

//	     &               status='old')
//	c chen reads a comment line
//	         read(20,*)
//	         read(20,*) ibinme,mapnum
        int ibinme_unknown;
        int mapnum_numOfLinesToBeRead;//this is the number of lines to be read from the beginning of fort21 file to
        //map_unknown (50 lines).
        Scanner scanner = new Scanner (data_from_fort21Reader_20.readLine());
        ibinme_unknown = scanner.nextInt();//this is 20
        mapnum_numOfLinesToBeRead = scanner.nextInt();//this is 50
        scanner.close();


        //----------------------------------------------------------------------------3 11.3.15



//		         if(ibinme.gt.ibin) then
//		         write(*,*) 'bin exceeds limit:',ibin,'!'
//		         stop
//      		 endif
        if(ibinme_unknown>ibin_unknown) throw new RuntimeException("bin exceeds limit: "+ibin_unknown+"!");


//		         do i=1,mapnum
//		         read(20,*) map(i)
//		         enddo

        for (int lineNumber=1;lineNumber<=mapnum_numOfLinesToBeRead;lineNumber++){
            scanner = new Scanner(data_from_fort21Reader_20.readLine());
            map_unknown[lineNumber]=scanner.nextInt();
            scanner.close();
        }


//
//		c tommer all ctmp are not used. they just enable skipping the beggining of the line with the names of the getAtoms.
//		 3001   read(20,*,end=3000) ctmp1,ctmp2,ctmp3,m,pp,i,j,k,l
//		c            if(pp.gt.9) pp=0.0
//		        if(i.eq.21.and.k.eq.21) then
//		c        hdpot(j,l)=pp
//		        else
//		c tommer fills the potential according to pp (pair potential). there is symmetry with respect to i,k because the order can be reversed.
//		        pot(i,j,k,l,m)=pp
//		        pot(k,l,i,j,m)=pp
//		        endif
//		        goto 3001


        //long startTimeforPotential = System.currentTimeMillis();
        class FortranArray5Dim_pot{
            float[][][][][] potential_pot=new float[ibin_unknown][20][15][20][15];//last index made first! 31.3.15
            protected float
            get_potential(int i1,int i2, int i3, int i4, int i5){
                return potential_pot[i5-1][i1-1][i2-1][i3-1][i4-1];//last index made first! 31.3.15
            }
            protected void
            set_potential(int i1,int i2, int i3, int i4, int i5, float value){
                potential_pot[i5-1][i1-1][i2-1][i3-1][i4-1]=value;//last index made first! 31.3.15
            }
        }

        FortranArray5Dim_pot pot_potential = new FortranArray5Dim_pot();

        String lineFromFortReader=data_from_fort21Reader_20.readLine();
        while(lineFromFortReader!=null){
            scanner = new Scanner(lineFromFortReader);
            int m1, i2, j3, k4, l5;
            float pairPotential_pp;
            scanner.next(); scanner.next(); scanner.next();
            m1=scanner.nextInt();
            pairPotential_pp=scanner.nextFloat();
            i2=scanner.nextInt();
            j3=scanner.nextInt();
            k4=scanner.nextInt();
            l5=scanner.nextInt();
            scanner.close();//done reading line
            if(!((i2==21)&&(k4==21))){//is this necessary? tommer 12.3.15
                pot_potential.set_potential(i2, j3, k4, l5, m1, pairPotential_pp);
                pot_potential.set_potential(k4, l5, i2, j3, m1, pairPotential_pp);//symmetry
            }

            lineFromFortReader=data_from_fort21Reader_20.readLine();
        }
        data_from_fort21Reader_20.close();



        //--------------------------------------------------------------------------------4 12.3.15



//		 3000   continue
//		        close(20)
//
//		              
//		c tommer this puts zeros in all of cnttheta's entries.
//		          do l1=1,20
//		          do m1=1,15
//		          do l2=1,20
//		          do m2=1,15
//		          do i=1,ibin
//		          do l=1,12
//		          do m=1,5
//		          cnttheta(l1,m1,l2,m2,i,l,m)=0.0
//		          enddo
//		          enddo
//		          enddo
//		          enddo
//		          enddo
//		          enddo
//		          enddo



        //long timeBeforeArray = System.currentTimeMillis();
        class FortranArray7Dim_cnttheta{
            float[][][][][][][] cnttheta=new float[20][15][20][15][ibin_unknown][5][12];//last 2 indices switched! 31.3.15
            protected float
            get_cnttheta(int i1,int i2, int i3, int i4, int i5, int i6, int i7){
                return cnttheta[i1-1][i2-1][i3-1][i4-1][i5-1][i7-1][i6-1];//last 2 indices switched! 31.3.15
            }
            protected void
            set_cnttheta(int i1,int i2, int i3, int i4, int i5, int i6, int i7, float value){
                cnttheta[i1-1][i2-1][i3-1][i4-1][i5-1][i7-1][i6-1]=value;//last 2 indices switched! 31.3.15
            }
        }
        FortranArray7Dim_cnttheta cnttheta_unknown=new FortranArray7Dim_cnttheta();

        //double[][][][][][][] cnttheta=new double[21][16][21][16][ibin+1][13][6];




//		        open(unit=21,file=base(1:ips)//'/fort.31_g72_noshift5_new',
//		     &               status='old')
        MeshiLineReader data_from_fort31Reader_21=
                new MeshiLineReader(Variables.base_goapalonePath+"\\fort.31_g72_noshift5_new");
//		         read(21,*)
        data_from_fort31Reader_21.readLine();//reads comment line
//		         read(21,*) ibinme,mapnum,ig 
        scanner = new Scanner(data_from_fort31Reader_21.readLine());
        ibinme_unknown=scanner.nextInt();
        mapnum_numOfLinesToBeRead=scanner.nextInt();
        int ig_s_parameter=scanner.nextInt();
        scanner.close();
//		         if(ibinme.gt.ibin) then
//		         write(*,*) 'bin exceeds limit:',ibin,'!'
//		         stop
//		         endif
        if(ibinme_unknown>ibin_unknown) throw new RuntimeException("bin exceeds limit: "+ibin_unknown+"!");
//		         do i=1,mapnum
//		         read(21,*) map(i)
//		c tommer this runs over all the previously read values in map() from the previous file. they might very well have the exact same value (should check!)
//		         enddo
        for (int i=1;i<=mapnum_numOfLinesToBeRead;i++){//this collects again values for map[].  (same values) 13.3.15
            scanner = new Scanner(data_from_fort31Reader_21.readLine());
            map_unknown[i]=scanner.nextInt();
            scanner.close();
        }


//
//
//		c tommer one in 6 lines is different. from there we get i,j,k,l.
//		 4001   read(21,*,end=4000) ctmp1,ctmp2,ctmp3,m,pp,i,j,k,l
//		c            if(pp.gt.9) pp=0.0
//		         do mm=1,5
//		c tommer from the other 5 lines we take m, mx and ppx(lx).
//		        read(21,*,end=4000) ctmp1,ctmp2,ctmp3,m,mx,
//		     &  (ppx(lx),lx=1,12),ppy,i,j,k,l
//		        if(i.eq.21.and.k.eq.21) then
//		        else
//		        do lx=1,12
//		        cnttheta(i,j,k,l,m,lx,mx)=ppx(lx)
//		        enddo 
//		        endif
//		         enddo 
//		       
//		c tommer again, symmetries (though these are a bit less straightforward).
//		        do lx=1,12
//		        cnttheta(k,l,i,j,m,lx,1)=cnttheta(i,j,k,l,m,lx,3)
//		        cnttheta(k,l,i,j,m,lx,2)=cnttheta(i,j,k,l,m,lx,4)
//		        cnttheta(k,l,i,j,m,lx,3)=cnttheta(i,j,k,l,m,lx,1)
//		        cnttheta(k,l,i,j,m,lx,4)=cnttheta(i,j,k,l,m,lx,2)
//		        cnttheta(k,l,i,j,m,lx,5)=cnttheta(i,j,k,l,m,lx,5)
//		        enddo 
//
//
//		        goto 4001
//		 4000   continue
//		        close(21)

        lineFromFortReader=data_from_fort31Reader_21.readLine();
        float[] pairPotentialX_ppx=new float[13];
        float pairPotentialY_ppy;
        int numberOfLineBeingRead=0;
        //long startTime = System.currentTimeMillis();
        int m,i,j,k,l, mx;
        int internalLineNumber, lx, placeInLine;
        float pairPotential_pp;
        data_from_fort31Reader_21.close();//for using File
        File fileFort31 = new File(Variables.base_goapalonePath+"\\fort.31_g72_noshift5_new");
        scanner = new Scanner(new FileReader(fileFort31));
        for (int lineToSkipNum=1;lineToSkipNum<=52;lineToSkipNum++) scanner.nextLine();
        String line = scanner.nextLine();
        String[] lineElements = new String[22];
        String cct_atomCode2letters;
        double[] xn2_unknown=new double[4], xd2_unknown= new double[4];
        double[] vt_normKtoJ= new double[4], vt2_normJtoK= new double[4];
        while(line!=null){
            boolean flag_noLine=false;
            numberOfLineBeingRead++;

            lineElements=line.split("\\s+");

            try{m=Integer.parseInt(lineElements[4]);}
            catch(Exception e){
                System.out.println("line count: "+numberOfLineBeingRead);
                for (int elementNum=0;elementNum<lineElements.length;elementNum++)
                    System.out.println("element number "+elementNum+" : "+lineElements[elementNum]);
                throw new RuntimeException(e);

            }
            pairPotential_pp=Float.parseFloat(lineElements[5]);
            try{i=Integer.parseInt(lineElements[6]);}
            catch(Exception e){
                System.out.println("line count: "+numberOfLineBeingRead);
                for (int elementNum=0;elementNum<lineElements.length;elementNum++)
                    System.out.println("element number "+elementNum+" : "+lineElements[elementNum]);
                throw new RuntimeException(e);

            }
            //i=Integer.parseInt(lineElements[6]);
            j=Integer.parseInt(lineElements[7]);
            k=Integer.parseInt(lineElements[8]);
            l=Integer.parseInt(lineElements[9]);

            if (scanner.hasNextLine()) line = scanner.nextLine();
            else {
                System.out.println("line: "+numberOfLineBeingRead );
                break;
            }
            for(internalLineNumber=1;(internalLineNumber<=5)&&(line!=null);internalLineNumber++){
                numberOfLineBeingRead++;
                lineElements=line.split("\\s+");
                mx =Integer.parseInt(lineElements[5]);
                for(placeInLine=1;placeInLine<=12;placeInLine++){
                    pairPotentialX_ppx[placeInLine]=Float.parseFloat(lineElements[5+placeInLine]);
                }
                pairPotentialY_ppy=Float.parseFloat(lineElements[18]);
                //INDICES have been changed (23.3.15). must be CHECKED!



                if(!((i==21)&&(k==21))){//is this necessary? tommer 13.3.15
                    for(lx=1;lx<=12;lx++){
                        cnttheta_unknown.set_cnttheta(i,j,k,l,m,lx,mx, pairPotentialX_ppx[lx]);
                    }
                }
                if (scanner.hasNextLine()) line = scanner.nextLine();
                else {
                    Utils.println("line: "+numberOfLineBeingRead+" , inside loop. inside index: "+internalLineNumber );
                    flag_noLine=true;
                    break;
                }
            }
            if (flag_noLine){
                if(internalLineNumber!=5) throw new RuntimeException("file ends with wrong line pattern (not 1+5)");
                break;
            }
            for(lx=1;lx<=12;lx++){//symmetries
                cnttheta_unknown.set_cnttheta(k,l,i,j,m,lx,1,cnttheta_unknown.get_cnttheta(i,j,k,l,m,lx,3));
                cnttheta_unknown.set_cnttheta(k,l,i,j,m,lx,2,cnttheta_unknown.get_cnttheta(i,j,k,l,m,lx,4));
                cnttheta_unknown.set_cnttheta(k,l,i,j,m,lx,3,cnttheta_unknown.get_cnttheta(i,j,k,l,m,lx,1));
                cnttheta_unknown.set_cnttheta(k,l,i,j,m,lx,4,cnttheta_unknown.get_cnttheta(i,j,k,l,m,lx,2));
                cnttheta_unknown.set_cnttheta(k,l,i,j,m,lx,5,cnttheta_unknown.get_cnttheta(i,j,k,l,m,lx,5));

            }



        }




        //--------------------------------------------------------------------------------5 13.3.15


//		          
//			open(unit=10,file=base(1:ips)//'/charge_inp.dat',
//		     &                     status='old')

        MeshiLineReader data_from_chargeReader_10=
                new MeshiLineReader(Variables.base_goapalonePath+"\\charge_inp.dat");
//		c tommer read comment line
//		        read(10,*) 
        data_from_chargeReader_10.readLine();
//		        id=1

        //declaring variables:
        final int ntyp_numberOfAminoAcids = 20;
        String[] resn_threeLetterResidueCodes = new String[ntyp_numberOfAminoAcids+1];
        int[] ianum_numbersOfAtomsInResidues= new int[21];
        String[][] cind_atomsByOrder = new String[ntyp_numberOfAminoAcids+1][21];
        for (int rowInd=1;rowInd<=ntyp_numberOfAminoAcids;rowInd++)
            for (int colInd=1;colInd<=20;colInd++)
                cind_atomsByOrder[rowInd][colInd]="XOP";
        float[][] chg_chargesByOrder = new float[ntyp_numberOfAminoAcids+1][16];//why 15 here and 20 above?

        String lineFromChargeReader;
        for(int id_residueIndex=1;data_from_chargeReader_10.ready();id_residueIndex++){
            //		c tommer read residue name and number of getAtoms in it (?)
            //		5        read (10,'(a3,1x,i2)',end=100) resnm,mm
            lineFromChargeReader=data_from_chargeReader_10.readLine();
            //		scanner = new Scanner(lineFromChargeReader);
            //		String residueCode_resm=scanner.next();
            //		int numOfAtomsInRes_mm=scanner.nextInt();
            lineElements = lineFromChargeReader.split("\\s+");
            String residueCode_resm=lineElements[0];
            int mm_numOfAtomsInRes=Integer.parseInt(lineElements[1]);




            //		        resn(id)=resnm
            resn_threeLetterResidueCodes[id_residueIndex]=residueCode_resm;


            //		        ianum(id)=mm
            ianum_numbersOfAtomsInResidues[id_residueIndex] = mm_numOfAtomsInRes;
            //		        do i=1,mm

            for(int i_atomNumInd = 1;i_atomNumInd<=mm_numOfAtomsInRes;i_atomNumInd++){
                //		c tommer read atom name and charge
                //		        read(10,*) atnm,xx
                lineFromChargeReader=data_from_chargeReader_10.readLine();
                lineElements = lineFromChargeReader.split("\\s+");

                String atnm_atomCode_tmp=lineFromChargeReader.substring(0,3);//!!!@@@
                //String atnm_atomCode_tmp=lineElements[0];//!!!@@@
                float xx_charge=Float.parseFloat(lineElements[1]);
                //		c tommer array of atomnames by residuenumber and atomnumber
                //		        cind(id,i)=atnm
                cind_atomsByOrder[id_residueIndex][i_atomNumInd]=atnm_atomCode_tmp;
                //		c tommer array of atomCharges by residuenumber and atomnumber
                //		        chg(id,i)=xx
                chg_chargesByOrder[id_residueIndex][i_atomNumInd]= xx_charge;
                //		        enddo
                //		        id=id+1
                //		        goto 5
            }
        }
//		100     continue
//		        close(10)
//
//
//

        //--------------------------------------------------------------------------------- 6 31.3.15
//		c chen Hooooooray the program starts
        //System.out.println("-------------programStart-----------\n\n\n");



        String filename;
        char chain;
        double ect_dFIRE, ect2_goapAG , ectTemp, ect2Temp;
        final int maxa_arraySize = 20001, maxres=1001;
        int[] resnum_firstLineOfResInPdbFileIndMinus1 = new int[maxa_arraySize];
        String whole_unknown;
        String resname_residueName; //char*4
        String[] whole1_linesFromFileArray = new String[maxa_arraySize];
        Double prob_errorCharacteristicNumber;
        String ares1_stringResidueNumberFromFile;
        int[][] ibk_idByResPlaceInFileAndAtomPlaceInRes = new int[maxres][16];
        int natom_unknown;
        double[] xp_xCoordinate= new double[maxa_arraySize], yp_yCoordinate= new double[maxa_arraySize], zp_zCoordinate= new double[maxa_arraySize+1];
        String atmname_atomCode;
        String[] rname_residueNamesByPlaceInFile = new String[maxa_arraySize];
        String[] aname_atomCodesFromFile = new String[maxa_arraySize+1];
        int[] restyp_residueCodeIndicesByPlaceInFile = new int[maxa_arraySize];
        int[] ihflg_flagForNotHydrogen = new int[maxa_arraySize];
        int[] ind1_resCodeIndexByI = new int[maxa_arraySize];

        int[] ind2_atomInResIndexByLine = new int[maxa_arraySize];
        float[] qq_unknown = new float[maxa_arraySize];
        double[][] xn_vectorList=new double[maxa_arraySize][4], xd_vectorList=new double[maxa_arraySize][4];
        double xxh_unknown[]=new double[5],yyh_unknown[]=new double[5],zzh_unknown[]=new double[5];//used before call for dihedral
        int[] ib0_unknown = new int[maxa_arraySize];
        int printCounter=1;


//		c A loop over all the model
        for(int ii_fileNumber=1;ii_fileNumber<=nfil_numberOfFiles;ii_fileNumber++){
//		         do ii=1,nfil
//		c            filename='/home/hzhou2/sparks/1010db/'//afil(ii)
//		            filename=afil(ii) THIS DOESN'T WORK RIGHT. I've changed it in the Fortran version.
            filename = Variables.base_goapalonePath+"/"+afil_filenames[ii_fileNumber];
//					
//		            if(chain.eq.'0') chain=' '
            chain = ' ';
//		c tommer ect = dFIRE? ect2=goapAG?
//		       ect=0.0
//		       ect2=0.0
            ect_dFIRE=0;
            ect2_goapAG=0;
//		c tommer this opens a pdb structure file
//		 1    open(11,file=filename,status='old') //this is seemingly NOT the start of a loop. 1.4.15

            scanner = new Scanner (new File(filename));
//		          ime=-1 
            int ime_redundantFlag=-1;
//		c tommer id is set to 1!
//		       id=1
            int id_lineInFileIndex=1;
//		       cx=0
//		       cy=0
//		       cz=0
//		        ires0=0
//		        ires=0
//		        ares0='     '
//		        rnum=0
//		        resnum(1)=0
            double cx_coordinateSum=0,cy_coordinateSum=0,cz_coordinateSum=0;
            int ires0_unknown=0, ires_residueInFileIndex=0, rnum_residuePlaceInFile=0;
            resnum_firstLineOfResInPdbFileIndMinus1[1]=0;
            String ares0_stringResidueInFileNumber = "     ";
//		c tommer it doean't seem like the variable "line" is used!!! there isn't even such a variable declared. 
//		c there is a character variable line1 which is declared but also not used. whatis the meaning of this code line?
//		c does this enable going to the next line in an open file?
//		 15    line = line + 1 //LOOP START!! 1.4.15 "line"????
            while(scanner.hasNextLine()){//LOOP to read from PDB file!
//		      
//		      read(11,'(a80)',err=10,end=10) whole1(id)
                //err=10 - jumps tp 10 if there is an error on read. 1.4.15
                whole1_linesFromFileArray[id_lineInFileIndex]=scanner.nextLine();
//		c      if(whole1(id)(22:22).ne.chain.and.ime.eq.1) goto 10 ////////THIS IS THE ONLY PLACE IME_REDUNDANTFLAG IS USED, and it's a comment.
//
//		c tommer to read the next line if doesn't start with 'ATOM'?
//		      if(whole1(id)(1:4).ne.'ATOM') goto 15
                if (!whole1_linesFromFileArray[id_lineInFileIndex].startsWith("ATOM"))
                    continue;
//
//		        
//		c      if(whole1(id)(1:4).ne.'ATOM'.or.whole1(id)(22:22).ne.chain)
//		c     &                     goto 15
//		         read(whole1(id)(18:21),'(a4)') resname 
                resname_residueName=whole1_linesFromFileArray[id_lineInFileIndex].substring(17, 21);
                //right indices?
//		             iusd=-1
                int iusd_flagForStartingWithResidueCode = -1;

//		          do i=1,20
//		          if(resn(i)(1:3).eq.resname(1:3)) iusd=1
//		          enddo
                for(int residueIndex=1;residueIndex<=20;residueIndex++){
                    if(resn_threeLetterResidueCodes[residueIndex].substring(0, 3).
                            equals(resname_residueName.substring(0, 3)))
                        iusd_flagForStartingWithResidueCode=1;
                }
//		c tommer read the next line if the line is not one of 20 amino acid types
//		          if(iusd.lt.0) goto 15
                if (iusd_flagForStartingWithResidueCode<0) continue;
//
//		      ime=1
                ime_redundantFlag=1;
//
//		      if(whole1(id)(17:17).ne.' ') 
//		     &  READ (whole1(id)(57:60),'(f4.2)') prob   //'(f4.2)' - 4 digits, 2 after decimal
                if (whole1_linesFromFileArray[id_lineInFileIndex].charAt(16)!=' '){
                    prob_errorCharacteristicNumber=
                            Double.parseDouble(whole1_linesFromFileArray[id_lineInFileIndex].substring(56, 60));

//		        
//		c tommer checks for reading next line. 1.d-4 - double precision?
//		      if(whole1(id)(17:17).ne.' '.and.prob.lt.0.5)   goto 15
                    if (prob_errorCharacteristicNumber<0.5) continue;



//		      if(whole1(id)(17:17).eq.'B'.and.abs(prob-0.5).lt.1.d-4) //1.d-4 -> double precision 1*10^-4 (?)
                    //		     &    goto 15
                    if ((whole1_linesFromFileArray[id_lineInFileIndex].charAt(16)=='B')&&
                            ((prob_errorCharacteristicNumber-0.5)<0.0001)) //what happens here when we have
                        // 0.50-0.5?
                        continue;
                }
//
//		      READ (whole1(id)(23:27),'(a5)') ares1
                ares1_stringResidueNumberFromFile = whole1_linesFromFileArray[id_lineInFileIndex].substring(22, 27);
//		c tommer if read character variable ares1 isn't '     ', fill column ires+1 in array ibk with -1.
//		c tommer ires starts from 0.
//		        if(ares1.ne.ares0) then
//		           ires=ires+1
//		           ares0=ares1
//		           do i=1,15
//		           ibk(ires,i)=-1
//		           enddo
//		        endif
                /*finalFeedback.println("id: "+id_lineInFileIndex+", ares1 :"+ares1_stringResidueNumberFromFile+":, ares0 :"+
                ares0_stringResidueInFileNumber+":\nires: "+ires_residueInFileIndex);*/
                if(!ares1_stringResidueNumberFromFile.equals(ares0_stringResidueInFileNumber)){
                    ires_residueInFileIndex=ires_residueInFileIndex+1; //ires was initiated with 0.
                    ares0_stringResidueInFileNumber=ares1_stringResidueNumberFromFile;
                    for (int ibkColIndex=1;ibkColIndex<=15;ibkColIndex++) {
                        ibk_idByResPlaceInFileAndAtomPlaceInRes[ires_residueInFileIndex][ibkColIndex] = -1;
                    }
                }// makes sure ares0 and ares1 are equal. 1.4.15
                //distinguishes the switch from one residue to the next (marked by a row of -1's)
//
//
//		c tommer id becomes atom number natom. id starts with 1.
//
//		      natom=id
                natom_unknown=id_lineInFileIndex; //starts with 1
//		c ires becomes residue number rnum.
//		      rnum=ires
                rnum_residuePlaceInFile= ires_residueInFileIndex;
//		  
//		c tommer value of id is saved in array resnum at (ires+1). why..? from earlier: resnum(1)=0.
//		        resnum(rnum+1)=id        

                resnum_firstLineOfResInPdbFileIndMinus1[rnum_residuePlaceInFile+1]=id_lineInFileIndex;
//		            
//		c           write(*,'(a80)') whole1(id)
//
//		      READ (whole1(id)(31:54),'(3F8.3)') Xp(id),Yp(id),Zp(id)
                lineElements=whole1_linesFromFileArray[id_lineInFileIndex].split("\\s+");
                xp_xCoordinate[id_lineInFileIndex]= Float.parseFloat(lineElements[5]);//coordinates from pdb file
                yp_yCoordinate[id_lineInFileIndex]= Float.parseFloat(lineElements[6]);
                zp_zCoordinate[id_lineInFileIndex]= Float.parseFloat(lineElements[7]);
//		c tommer coordinates/potential? need to check pdb format.
//		       cx=cx+xp(id)
//		       cy=cy+yp(id)
//		       cz=cz+zp(id)
                cx_coordinateSum=cx_coordinateSum+xp_xCoordinate[id_lineInFileIndex];
                cy_coordinateSum=cy_coordinateSum+yp_yCoordinate[id_lineInFileIndex];
                cz_coordinateSum=cz_coordinateSum+zp_zCoordinate[id_lineInFileIndex];
//
//
//
//
//		         read(whole1(id)(14:17),'(a4)') atmname 
                atmname_atomCode = whole1_linesFromFileArray[id_lineInFileIndex].substring(13, 17);


//		         if(atmname(1:2).eq.'OT') atmname='O   '
                if (atmname_atomCode.substring(0, 2).equals("OT")) atmname_atomCode="O   ";
//		         if(atmname(1:3).eq.'OXT') atmname='O   '   
                if (atmname_atomCode.substring(0, 3).equals("OXT")) atmname_atomCode="O   ";
//		         if(resname(1:3).eq.'ILE'.and.atmname(1:3).eq.'CD1') 
//		     &           atmname='CD  '   
                if (resname_residueName.substring(0, 3).equals("ILE")&&
                        atmname_atomCode.substring(0, 3).equals("CD1"))
                    atmname_atomCode="CD  ";

//		         rname(rnum)=resname
                rname_residueNamesByPlaceInFile[rnum_residuePlaceInFile]=resname_residueName;

//		         aname(id)=atmname
                aname_atomCodesFromFile[id_lineInFileIndex]= atmname_atomCode;
//		          restyp(rnum)=0
                restyp_residueCodeIndicesByPlaceInFile[rnum_residuePlaceInFile]=0; //why??
//		      
//		c tommer i=id. why is this extra index needed???
//				 i=id

                int i_idDUPLICATEmaybe= id_lineInFileIndex;
//
//
//		c tommer flag for "is not hydrogen"
//		            ihflg(id)=-1
                ihflg_flagForNotHydrogen[id_lineInFileIndex]=-1;
//		            ind2(i)=-1
                ind2_atomInResIndexByLine[i_idDUPLICATEmaybe]=-1;

//		c         if(atmname(1:1).ne.'H'.and.atmname(1:1).ne.'A') then

//		         if(atmname(1:1).ne.'H')then
//		            ihflg(id)=1
                if (atmname_atomCode.charAt(0)!='H'){ //BIG IF! 27.4.15
                    ihflg_flagForNotHydrogen[id_lineInFileIndex]=1;

                    //		         iok=-1
                    int iok_flagIdentifiedAtom=-1;
                    //		         do j1=1,20
                    //		         if(resn(j1).eq.resname(1:3)) then
                    //		c tommer to get from rnum to residue type it is necessary to use restyp(rnum).
                    //		            restyp(rnum)=j1
                    //		            do j2=1,15
                    //		            if(cind(j1,j2).eq.atmname(1:3)) then
                    //		c            qq(i)=chg(j1,j2)
                    //		            ind1(i)=j1
                    //		            ind2(i)=j2
                    //		c tommer ibk seems to contain the atom id-s ordered by rnum and place within residue.
                    //		            ibk(rnum,j2)=id
                    //		c check that non-hydrogen atom is identified ok.
                    //		            iok=1
                    //		            endif
                    //		            enddo
                    //		         endif
                    //		         enddo
                    for(int j1_residueCodeIndex=1;j1_residueCodeIndex<=20;j1_residueCodeIndex++){
                        if (resn_threeLetterResidueCodes[j1_residueCodeIndex].equals
                                (resname_residueName.substring(0, 3))){
                            restyp_residueCodeIndicesByPlaceInFile[rnum_residuePlaceInFile]=j1_residueCodeIndex;
                            for (int j2_atomInResIndex=1;j2_atomInResIndex<=15;j2_atomInResIndex++){
                                //System.out.println("residueCodeIndex: "+residueCodeIndex+", atomInResInd: "+atomInResIndex);
                                if(cind_atomsByOrder[j1_residueCodeIndex][j2_atomInResIndex]==null) {

                                    scanner.close();
                                    throw new RuntimeException("cind entry is null, cannot compare");
                                }
                                //THIS IS TO PREVENT NULL EXCEPTION. this is NOT how it's done in the original code. 27.4.15








                                if (cind_atomsByOrder[j1_residueCodeIndex][j2_atomInResIndex].equals
                                        (atmname_atomCode.substring(0, 3))){
                                    ind1_resCodeIndexByI[i_idDUPLICATEmaybe]=j1_residueCodeIndex;
                                    ind2_atomInResIndexByLine[i_idDUPLICATEmaybe]=j2_atomInResIndex;
                                    ibk_idByResPlaceInFileAndAtomPlaceInRes[rnum_residuePlaceInFile][j2_atomInResIndex]=id_lineInFileIndex; //here i==id...?
                                    //atom in residue is always the same for every residue of the same type. therefore, this array
                                    //gives the line in file by the number of the residue in the file and the number of the atom in that residue.
                                    iok_flagIdentifiedAtom=1;


                                }
                            }
                        }
                    }
                    //feedbackWriter.close();

                    //
                    //		c tommer atom not found. why chose these indices? (10,1)? this looks like a bug. wouldn't disturb the program on normal input...
                    //		         if(iok.lt.0) then
                    //		c         write(*,*) afil(ii),ires,resname(1:3),atmname(1:3),
                    //		c     &            'not found !'
                    //		            ind1(i)=10
                    //		            ind2(i)=1
                    //		            ibk(rnum,1)=id
                    //		c         stop
                    //		         endif

                    if (iok_flagIdentifiedAtom<0){

                        ind1_resCodeIndexByI[i_idDUPLICATEmaybe]=10;// this is 'GLY'?
                        ind2_atomInResIndexByLine[i_idDUPLICATEmaybe]=1;//this is 'N' for Nitrogen??
                        //why GLY Nitrogen?
                        ibk_idByResPlaceInFileAndAtomPlaceInRes[rnum_residuePlaceInFile][1]=id_lineInFileIndex;
                        //why insert the line number in the place of the first atom?
                    }
                    //
                }
//		         endif 
//
//
//		      natom=id
                natom_unknown=id_lineInFileIndex; //isn't this already the case?? [aprox. 120 lines before...]
//		c tommer end of loop to read pdb file
//		      id=id+1
                id_lineInFileIndex++;
//		      goto 15

            } //end of 15 loop
//
//
//		 10   continue
//
//		           close(11)
            scanner.close();
//

            /*finalFeedback.println("resnum:");
            for (int i_runner=1; i_runner<20;i_runner++)
                finalFeedback.println(i_runner+": "+resnum_firstLineOfResInPdbFileIndMinus1[i_runner]);*/


            int ibb_unknown=0;
//		          do k1=2,rnum
            //new inner loop. goes over residues by order.
            //variables are still read from the pdb file! 4.5.15
            for (int k1_residueInd=2;k1_residueInd<=rnum_residuePlaceInFile;k1_residueInd++){
//		         if((resnum(k1+1)-resnum(k1)).lt.ianum(restyp(k1))) goto 404

                if ((resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd+1]-resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd])<
                        ianum_numbersOfAtomsInResidues[restyp_residueCodeIndicesByPlaceInFile[k1_residueInd]])
                    continue;//checks that there is the right number of getAtoms in the residue.




//		           do k=resnum(k1)+1,resnum(k1+1) 
				/*System.out.println("resnum[1]: "+resnum_firstLineOfResInPdbFileIndMinus1[1]+
						"\nresnum[2]: "+resnum_firstLineOfResInPdbFileIndMinus1[2]+
						"\nresnum[3]: "+resnum_firstLineOfResInPdbFileIndMinus1[3]);*/
                for (int k_lineInFileIndex_inResidue=resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd]+1;
                     k_lineInFileIndex_inResidue<=resnum_firstLineOfResInPdbFileIndMinus1[k1_residueInd+1];
                     k_lineInFileIndex_inResidue++){
//		           if(ind2(k).gt.0) then



                    if(ind2_atomInResIndexByLine[k_lineInFileIndex_inResidue]>0)	//ind2 was given a value
//		           call caldplane(xp,yp,zp,k1,k,ind1(k),ind2(k),ibk,xn1,xd1,ibb)



                        //***************************************now we call CALDPLANE ***********************************************8

//                    if(k_lineInFileIndex_inResidue==10){
//                        finalFeedback.println("\nxp: "+xp_xCoordinate+", "+yp_yCoordinate+", "+zp_zCoordinate+
//                                "\nk1_residueInd: "+k1_residueInd+", k_lineinFilIndex_inResidue: "+k_lineInFileIndex_inResidue+
//                        "\n ind1[k]: "+ind1_resCodeIndexByI[k_lineInFileIndex_inResidue]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInFileIndex_inResidue]+"\n\n");
//                    }

                        ib0_unknown[k_lineInFileIndex_inResidue]=caldplane_unknown_clean(xp_xCoordinate,yp_yCoordinate,
                                zp_zCoordinate,
                                k1_residueInd,k_lineInFileIndex_inResidue,ind1_resCodeIndexByI[k_lineInFileIndex_inResidue],
                                ind2_atomInResIndexByLine[k_lineInFileIndex_inResidue],ibk_idByResPlaceInFileAndAtomPlaceInRes);/*,
								Variables.xn1_VectorInput_unknown,Variables.xd1_VectorInput_unknown*///finalFeedback);//ibb is given a value WITHIN THE SUBROUTINE. apparently, unlike in a java function,
                    //this value is preserved when leaving the subroutine. Tommer 18.4.15

//		           ib0(k)=ibb  //this is implemented above


                    //--------tentative explanations about variables in the function------------
					/*the function gets as input:
					 * 1. coordinate vectors (xp[], yp[], zp[])
					 * 2. place of residue in pdb file [first residue, second residue etc.] (k1)
					 * 3. line-number  of first line in residue (k)
					 * 4. indices of atom and residue code (ind2 and ind1)
					 * 
					 * 
					 * (internal variables) or RETURNED VARIABLES:
					 * 
					 * 
					 * (1.) v2: i3-i1, vector from next atom to atom or from previous previous atom to atom
					 * (2. additional notation-not a variable) v4: vector from previous atom atom to atom
					 * 3. xn1: v2+v4 or v4
					 * 4. xd1: (v2Xv4)X(v2+v4) or (xn1Xv2)X xn1
					 * 
					 * 5. ibb: the place in the pdb file 
					 * 
					 * 
					 * ADDITIONAL NOTEWORTHY RELATIONS:
					 * if the above descriptions are accurate, it follows that always:
					 * - xn1 and xd1 are orthogonal
					 * - xd1 and (v2+v4) are orthogonal
					 */




                    //if (1==1) throw new RuntimeException("this is the end of the checkup");
                    //----------------------------------------------------------------------------------------------------------------------
                    //------------------------END OF CHECKUP 27.5.15-------------------------------------------//


//			           do l=1,3
                    for(int axisInd=1;axisInd<=3;axisInd++){
//			           xn(k,l)=xn1(l)
//			           xd(k,l)=xd1(l)
                        xn_vectorList[k_lineInFileIndex_inResidue][axisInd]=Variables.xn1_Vz_unnormalized[axisInd];
                        xd_vectorList[k_lineInFileIndex_inResidue][axisInd]=Variables.xd1_Vx_unnormalized[axisInd];

//			           enddo 
                    }// create a list of all acquaired xn1, xd1's.
//			           endif
//			           
                }
//				   enddo
            }
//		404	continue
//		          enddo
//
//            finalFeedback.println("xn list:");
//            for (int k_ind=1;k_ind<=100;k_ind++) {
//                finalFeedback.println("xn["+k_ind+"]: " + xn_vectorList[k_ind][1]+", "+xn_vectorList[k_ind][2]+", "+xn_vectorList[k_ind][3]);
//            }
//            finalFeedback.println("\n\n");
//
//		         do k1=2,rnum-1
            for (int k1_resInd=2; k1_resInd<=rnum_residuePlaceInFile-1;k1_resInd++){
                //loop over residues, from second in pdb file(i think) to 1 before rnum'th.
//		         if((resnum(k1+1)-resnum(k1)).lt.ianum(restyp(k1))) goto 405
                if((resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd+1]-
                        resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd])<ianum_numbersOfAtomsInResidues
                        [restyp_residueCodeIndicesByPlaceInFile[k1_resInd]]) continue;
                //check that there's the right number of getAtoms in the residue (bug control?)
                //

//		           do k=resnum(k1)+1,resnum(k1+1)   
                for (int k_lineInRes=resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd]+1;k_lineInRes<=
                        resnum_firstLineOfResInPdbFileIndMinus1[k1_resInd+1];k_lineInRes++){
                    //loop over the lines of the k1'th residue

//		           cct=aname(k)(1:2)
                    cct_atomCode2letters=aname_atomCodesFromFile[k_lineInRes].substring(0, 2);
//		           if(ind2(k).gt.0) then
                    if(ind2_atomInResIndexByLine[k_lineInRes]>0){//file not over
                        //if 1.
//		           do l=1,3
                        for(int l_dimInd=1;l_dimInd<=3;l_dimInd++){
//		           xn1(l)=xn(k,l)
//		           xd1(l)=xd(k,l)
                            Variables.xn1_Vz_unnormalized[l_dimInd]=xn_vectorList[k_lineInRes]
                                    [l_dimInd];
                            Variables.xd1_Vx_unnormalized[l_dimInd]=xd_vectorList[k_lineInRes]
                                    [l_dimInd];
                        }//regain our local xn1, xd1
//		           enddo


//
//
//		             do j1=k1+1,rnum
                        for (int j1_resCurrentToLast=k1_resInd+1;j1_resCurrentToLast<=//j1key
                                rnum_residuePlaceInFile;j1_resCurrentToLast++)
                        {

//		         if((resnum(j1+1)-resnum(j1)).lt.ianum(restyp(j1))) goto 406

                            if((resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast+1]-
                                    resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast])<
                                    ianum_numbersOfAtomsInResidues
                                            [restyp_residueCodeIndicesByPlaceInFile[j1_resCurrentToLast]])
                                //	break;//zaq4 THIS CAN BREAK THE LOOP AND MAKE TROUBLE. Indeed. Chen replaced by continue;
                                continue;
//		             do j=resnum(j1)+1,resnum(j1+1)
                            //if (1==1) throw new RuntimeException("j1: "+j1_resCurrentToLast+", resnum[j1]: "+
                            //resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast]);
                            for(int j_lineInRes=resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast]+1;
                                j_lineInRes<=resnum_firstLineOfResInPdbFileIndMinus1[j1_resCurrentToLast+1];
                                j_lineInRes++)
                            {
//		             cct=aname(j)(1:2)
                                cct_atomCode2letters=
                                        aname_atomCodesFromFile[j_lineInRes].substring(0, 2);
//		           if(ind2(j).gt.0.) then






                                if(ind2_atomInResIndexByLine[j_lineInRes]>0){//file not over

                                    //if 2.

//		           do l=1,3
                                    for(int l_dim=1;l_dim<=3;l_dim++){
//		           xn2(l)=xn(j,l)
//		           xd2(l)=xd(j,l)
                                        xn2_unknown[l_dim]=xn_vectorList[j_lineInRes][l_dim];
                                        xd2_unknown[l_dim]=xd_vectorList[j_lineInRes][l_dim];
//		           enddo
                                    }


//
//		              rd=sqrt((xp(k)-xp(j))**2+(yp(k)-yp(j))**2+
//		     &                (zp(k)-zp(j))**2)
                                    double rd_distanceKJ = Math.sqrt(
                                            (xp_xCoordinate[k_lineInRes]-xp_xCoordinate[j_lineInRes])*(xp_xCoordinate[k_lineInRes]-xp_xCoordinate[j_lineInRes])+
                                                    (yp_yCoordinate[k_lineInRes]-yp_yCoordinate[j_lineInRes])*(yp_yCoordinate[k_lineInRes]-yp_yCoordinate[j_lineInRes])+
                                                    (zp_zCoordinate[k_lineInRes]-zp_zCoordinate[j_lineInRes])*(zp_zCoordinate[k_lineInRes]-zp_zCoordinate[j_lineInRes]));

//
//
//
//		c             if(jj.gt.ibin) jj=ibin 
//
//		             jj=map(int(rd*2))
                                    int jj_unknown=map_unknown[(int)(rd_distanceKJ*2)];//DOES CASTING WORK THE SAME WAY? 27.5
//
//		             if(jj.le.ibin.and.jj.gt.0.1.and.rd.lt.30) then 







                                    if ((jj_unknown<=ibin_unknown)&&(jj_unknown>0.1)&&(rd_distanceKJ<30)){
                                        //unknown boubdary check...
                                        //if 3.



//
//		             ee=pot(ind1(k),ind2(k),ind1(j),ind2(j),jj)
//		             ect=ect+ee
                                        double ee_KJjjPotential=pot_potential.get_potential(
                                                ind1_resCodeIndexByI[k_lineInRes],
                                                ind2_atomInResIndexByLine[k_lineInRes],
                                                ind1_resCodeIndexByI[j_lineInRes],
                                                ind2_atomInResIndexByLine[j_lineInRes],
                                                jj_unknown);
                                        ectTemp=ect_dFIRE;//checkup 8.6.15
                                        ect_dFIRE=ect_dFIRE+ee_KJjjPotential;



//
//		             if((j1-k1).ge.ig) then




                                        if((j1_resCurrentToLast-k1_resInd)>=ig_s_parameter){//only proceed
                                            //if distance between residues is higher or equal to s parameter (7)
                                            //for closer residues angular part is not random
                                            // if 4.
//		             xdd=1./rd
                                            double xdd_distReciprocal= 1./rd_distanceKJ;
//		             vt(1)=(xp(j)-xp(k))*xdd   ! 1->2
//		             vt(2)=(yp(j)-yp(k))*xdd   ! 1->2
//		             vt(3)=(zp(j)-zp(k))*xdd   ! 1->2
                                            vt_normKtoJ[1] = (xp_xCoordinate[j_lineInRes]-
                                                    xp_xCoordinate[k_lineInRes])*xdd_distReciprocal;
                                            vt_normKtoJ[2]=(yp_yCoordinate[j_lineInRes]-
                                                    yp_yCoordinate[k_lineInRes])*xdd_distReciprocal;
                                            vt_normKtoJ[3]=(zp_zCoordinate[j_lineInRes]-
                                                    zp_zCoordinate[k_lineInRes])*xdd_distReciprocal;

//
//		             vt2(1)=-vt(1)
//		             vt2(2)=-vt(2)
//		             vt2(3)=-vt(3)
                                            vt2_normJtoK[1]=-vt_normKtoJ[1];
                                            vt2_normJtoK[2]=-vt_normKtoJ[2];
                                            vt2_normJtoK[3]=-vt_normKtoJ[3];




                                            //feedbackWriter.close();
                                            //if (1==1) throw new RuntimeException("this is the end of the checkup 27.5");
//

                                            //-----------------------------------------------------------27.5.15-------------------------------------------------------------------

//		             call calang2(xn1,vt,ang1,cs1)
                                            double  pointerToAngAndCos1_thetaA_cosThetaA[]=new double[3], pointerToAngAndCos2_thetaB_cosThetaB[]=new double[3];
                                            int mm1_cosThetaA_xn1_vt_bin, mm2_psyA_xd1_vt_bin, mm3_cosThetaB_xn2_vt2_bin, mm4_psyB_xd2_vt2_bin;

                                            calang2_calculateAngleInDeg(Variables.xn1_Vz_unnormalized,vt_normKtoJ, pointerToAngAndCos1_thetaA_cosThetaA);


//		             mm1=int((cs1+1.001)*6.0)+1
                                            mm1_cosThetaA_xn1_vt_bin= (int)((pointerToAngAndCos1_thetaA_cosThetaA[2]+1.001)*6.0)+1;
//		             if(mm1.gt.12) mm1=12
                                            if(mm1_cosThetaA_xn1_vt_bin>12) mm1_cosThetaA_xn1_vt_bin=12;
//		             call calphi(xn1,xd1,vt,phi1)
                                            double phi1_signedAngleBet_xd1Andvt_psyA=
                                                    calphi_signedAngleBet_xh1ANDxd_PSY(Variables.xn1_Vz_unnormalized,
                                                            Variables.xd1_Vx_unnormalized,vt_normKtoJ);

//		             mm2=int((phi1+180.001)/30.0)+1
                                            mm2_psyA_xd1_vt_bin= (int)((phi1_signedAngleBet_xd1Andvt_psyA+
                                                    180.001)/30.0)+1;
//		             if(mm2.gt.12) mm2=12
                                            if(mm2_psyA_xd1_vt_bin>12) mm2_psyA_xd1_vt_bin=12;
//----------------------------------------------------------------1.6.15 checkup----------------------------------------




                                            //feedbackWriter.close();
                                            //if (1==1) throw new RuntimeException("this is the end of the checkup 1.6");




//		             call calang2(xn2,vt2,ang2,cs2)
//                                           if(printCounter==1){
//                                                finalFeedback.println("\nxn2: "+xn2_unknown[1]+", "+xn2_unknown[2]+", "+xn2_unknown[3]+"\nvt2: "+
//                                                        vt2_normJtoK[1]+", "+vt2_normJtoK[2]+", "+vt2_normJtoK[3]+"\n");
//                                            }
                                            calang2_calculateAngleInDeg(
                                                    xn2_unknown, vt2_normJtoK, pointerToAngAndCos2_thetaB_cosThetaB);




//		             mm3=int((cs2+1.001)*6.0)+1
                                            mm3_cosThetaB_xn2_vt2_bin=(int)((pointerToAngAndCos2_thetaB_cosThetaB[2]+1.001)*6.0)+1;
//		             if(mm3.gt.12) mm3=12
                                            if(mm3_cosThetaB_xn2_vt2_bin>12) mm3_cosThetaB_xn2_vt2_bin=12;
//		             call calphi(xn2,xd2,vt2,phi2)
                                            double phi2_signedAngleBet_xd2Andvt2_psyB=calphi_signedAngleBet_xh1ANDxd_PSY(
                                                    xn2_unknown, xd2_unknown, vt2_normJtoK);
//		             mm4=int((phi2+180.001)/30.0)+1
                                            mm4_psyB_xd2_vt2_bin=(int)((phi2_signedAngleBet_xd2Andvt2_psyB+180.001)/30.0)+1;

//		             if(mm4.gt.12) mm4=12
                                            if(mm4_psyB_xd2_vt2_bin>12) mm4_psyB_xd2_vt2_bin=12;
//
//--------------------------------------------------1.6 checkup---------------------------------------




//			   		xxh(1)=xp(j)+xn2(1)
//		             yyh(1)=yp(j)+xn2(2)
//		             zzh(1)=zp(j)+xn2(3)
//		             xxh(2)=xp(j)
//		             yyh(2)=yp(j)
//		             zzh(2)=zp(j)
//		             xxh(3)=xp(k)
//		             yyh(3)=yp(k)
//		             zzh(3)=zp(k)
//		             xxh(4)=xp(k)+xn1(1)
//		             yyh(4)=yp(k)+xn1(2)
//		             zzh(4)=zp(k)+xn1(3)
                                            xxh_unknown[1]=xp_xCoordinate[j_lineInRes]+xn2_unknown[1];
                                            yyh_unknown[1]=yp_yCoordinate[j_lineInRes]+xn2_unknown[2];
                                            zzh_unknown[1]=zp_zCoordinate[j_lineInRes]+xn2_unknown[3];
                                            xxh_unknown[2]=xp_xCoordinate[j_lineInRes];
                                            yyh_unknown[2]=yp_yCoordinate[j_lineInRes];
                                            zzh_unknown[2]=zp_zCoordinate[j_lineInRes];
                                            xxh_unknown[3]=xp_xCoordinate[k_lineInRes];
                                            yyh_unknown[3]=yp_yCoordinate[k_lineInRes];
                                            zzh_unknown[3]=zp_zCoordinate[k_lineInRes];
                                            xxh_unknown[4]=xp_xCoordinate[k_lineInRes]+Variables.xn1_Vz_unnormalized[1];
                                            yyh_unknown[4]=yp_yCoordinate[k_lineInRes]+Variables.xn1_Vz_unnormalized[2];
                                            zzh_unknown[4]=zp_zCoordinate[k_lineInRes]+Variables.xn1_Vz_unnormalized[3];
//



//		            call dihedral(xxh,yyh,zzh,dih)
                                            double dih_angleBetweenNormsToConsecutivePlanes_xhi=
                                                    dihedral_calc_Xhi(xxh_unknown,yyh_unknown,zzh_unknown);

											
											
											
											
											
											
											
											/*if(printCounter==11172){//checkup 30.6.15
												System.out.println("checkup 30.6\n"+"dih: "+dih_angleBetweenNormsToConsecutivePlanes);
												System.out.println("xxh: "+xxh_unknown[1]+", "+xxh_unknown[2]+", "+xxh_unknown[3]+", "+xxh_unknown[4]);
												System.out.println("yyh: "+yyh_unknown[1]+", "+yyh_unknown[2]+", "+yyh_unknown[3]+", "+yyh_unknown[4]);
												System.out.println("zzh: "+zzh_unknown[1]+", "+zzh_unknown[2]+", "+zzh_unknown[3]+", "+zzh_unknown[4]);
												System.out.println("\n");
											}*/



//		             mm5=min(int((dih+180.001)/30.0)+1,12)
                                            int mm5_dih_xhi_bin=Math.min(
                                                    (int)((dih_angleBetweenNormsToConsecutivePlanes_xhi+180.001)/30.0)+1, 12);
//
//
//
//
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm1,1)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm2,2)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm3,3)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm4,4)
//		          ect2=ect2+cnttheta(ind1(k),ind2(k),ind1(j),ind2(j),jj,mm5,5)
                                            ect2Temp=ect2_goapAG;

                                            if(printCounter==1) {
                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm1_cosThetaA_xn1_vt_bin, 1);
//                                                finalFeedback.println("\nfirst addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp));

                                                ect2Temp=ect2_goapAG;
                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm2_psyA_xd1_vt_bin, 2);
//                                                finalFeedback.println("\nsecond addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp));

                                                ect2Temp=ect2_goapAG;
                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm3_cosThetaB_xn2_vt2_bin, 3);
//                                                finalFeedback.println("\nthird addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                "\nk: "+k_lineInRes+", j: "+j_lineInRes+
//                                                "\nind1[k]: "+ind1_resCodeIndexByI[k_lineInRes]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInRes]+
//                                                "\nind1[j]: "+ind1_resCodeIndexByI[j_lineInRes]+", ind2[j]: "+ind2_atomInResIndexByLine[j_lineInRes]+
//                                                "\njj: "+jj_unknown+", mm3: "+mm3_cosThetaB_xn2_vt2_bin);

                                                ect2Temp=ect2_goapAG;
                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm4_psyB_xd2_vt2_bin, 4);
//                                                finalFeedback.println("\nfourth addition"+
//													"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                "\nk: "+k_lineInRes+", j: "+j_lineInRes+
//                                                "\nind1[k]: "+ind1_resCodeIndexByI[k_lineInRes]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInRes]+
//                                                "\nind1[j]: "+ind1_resCodeIndexByI[j_lineInRes]+", ind2[j]: "+ind2_atomInResIndexByLine[j_lineInRes]+
//                                                "\njj: "+jj_unknown+", mm4: "+mm4_psyB_xd2_vt2_bin);
                                                ect2Temp=ect2_goapAG;
                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm5_dih_xhi_bin, 5);
//                                                finalFeedback.println("\nfifth addition"+
//												"\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                "\nk: "+k_lineInRes+", j: "+j_lineInRes+
//                                                "\nind1[k]: "+ind1_resCodeIndexByI[k_lineInRes]+", ind2[k]: "+ind2_atomInResIndexByLine[k_lineInRes]+
//                                                "\nind1[j]: "+ind1_resCodeIndexByI[j_lineInRes]+", ind2[j]: "+ind2_atomInResIndexByLine[j_lineInRes]+
//                                                "\njj: "+jj_unknown+", mm5: "+mm5_dih_xhi_bin);
                                            }
                                            else{
                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm1_cosThetaA_xn1_vt_bin, 1);

                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm2_psyA_xd1_vt_bin, 2);

                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm3_cosThetaB_xn2_vt2_bin, 3);

                                                ect2_goapAG = ect2_goapAG + cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm4_psyB_xd2_vt2_bin, 4);

                                                ect2_goapAG = ect2_goapAG +cnttheta_unknown.get_cnttheta(
                                                        ind1_resCodeIndexByI[k_lineInRes], ind2_atomInResIndexByLine[k_lineInRes],
                                                        ind1_resCodeIndexByI[j_lineInRes], ind2_atomInResIndexByLine[j_lineInRes],
                                                        jj_unknown, mm5_dih_xhi_bin, 5);

//                                                 finalFeedback.println("\nect2: "+ect2_goapAG+
//													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
//                                                    "\ncount: "+printCounter);
                                            }


											/*finalFeedback.println("\nect: "+ect_dFIRE+
													"\nect increment: "+(ect_dFIRE-ectTemp)+
													"\nect2: "+ect2_goapAG+
													"\nect2 increment: "+(ect2_goapAG-ect2Temp)+
													"\ncount: "+printCounter);*/




                                            printCounter++;




                                            //feedbackWriter.close();
                                            //if (1==1) throw new RuntimeException("this is the end of the checkup 3.6");


//		             endif ! .ge.ig
                                        }//if (j1-k1)> ig
//
//		             endif
                                    }// boundary check for jj and rd
//
//		            endif
                                }
//		1234              continue
//		            enddo


                            }
                            //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@zaq3@@@@@@@@@@@
                            //this is the last place where the fortran program
                            //can be exited without error. 6.5.15


//		406	continue
//		            enddo
                        }
//
//		           endif
                    }
//		          enddo
                }
//
//
//		405	continue
//		         enddo
            }
//
//		       write(*,555) ii,' ',afil(ii),ect+ect2,ect,ect2
//		403	continue
//		      enddo
//		c tommer this is the end of the main loop??
//
//		555	format(1x,i5,a1,1x,a30,1x,f12.2,1x,2f12.2)
//		     

            System.out.println("\n\n\nOutput, 30.6.15:\nii: "+ii_fileNumber+", afil[ii]: "+afil_filenames[ii_fileNumber]+
                    "\ngoap score: "+(ect_dFIRE+ect2_goapAG)+"\ndFIRE score: "+ect_dFIRE+"\ngoapAG score: "+ect2_goapAG);//checkup
            //finalFeedback.println("\n\n\nOutput, 6.6.15:\nii: "+ii_fileNumber+", afil[ii]: "+afil_filenames[ii_fileNumber]+
            //	"\ngoap score: "+(ect_dFIRE+ect2_goapAG)+"\ndFIRE score: "+ect_dFIRE+"\ngoapAG score: "+ect2_goapAG);//checkup 8.6.15



        } //end of model loop?-----------------------------------------------------------------------------------------

//        stop
//      end		

//		finalFeedback.close();

    }
    //end of main



    private static class Variables{//static variables for soubroutines
        static int ips_goapPathLength;//this is from "common"
        static String base_goapalonePath;//this is from "common"
				/*protected void setIps(int ips){
					ips_goapPathLength=ips;
				}
				protected void setBase(String base){
					base_goapalonePath=base;
				}
				protected int getIps() {
					return ips_goapPathLength;
				}
				protected String getBase(){
					return base_goapalonePath;
				}//end of "common" variables*/

        //variables for subroutine caldplane:
        static boolean init_flagFirstTimeIncaldplane=true;
        static int[][][] sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine = new int[21][4][16];
        static double[] xn1_Vz_unnormalized=new double[4], xd1_Vx_unnormalized=new double[4];

    }



    private int caldplane_unknown_cleaner(double xp_unknown[], double yp_unknown[], double[] zp_unknown,
                                          int k1_residueInd, int k_lineInFileIndex_inResidue, int id1_resCodeInd, int id2_atomInResInd,
                                          int[][] ibk_idByResPlaceInFileAndAtomPlaceInRes) {
        int ibb_previousAtomLine; //this will be returned.


        //double[][][] sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine = new double[21][5][16];
        //this means it will be created anew each time the function is called!
        //this variable is not used in the function!

        int l_lineLengthInGeometryFile;
        double[] v1_Vy_unnormalized= new double[4], v2_r13_unnormalized=new double[4];//, v3_unknown=new double[4];


        int i1_kReplica_lineInFileIndex = k_lineInFileIndex_inResidue, i2_previousAtomLineA1=0, i3_nextAtomLineA2_orPrevPrev=0;
//			         n0=-1
        int n0_flagForNextFound=-1;
//			         if(id2.eq.1) then ! N
        if(id2_atomInResInd==1){//Nitrogen //FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1-1,3)        ! k1-1's C
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd-1][3];
//			           i3=ibk(k1,2)          ! k1's CA
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][2];
//			           n0=1
            n0_flagForNextFound=1;
        }

//			         elseif(id2.eq.2) then !  CA
        if(id2_atomInResInd==2){ //FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1,1)         !k1  N
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][1];
//			           i3=ibk(k1,3)         !k1  C
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][3];
//			           n0=1
            n0_flagForNextFound=1;
        }

//			         elseif(id2.eq.3) then !  C
        if(id2_atomInResInd==3){//FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1,2)         !k1  CA
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][2];
//			           i3=ibk(k1,4)         !k1  O
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][4];
//			           n0=1
            n0_flagForNextFound=1;
        }

//			        elseif(id2.eq.4) then !  O
        if(id2_atomInResInd==4){//FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1,3)         !k1  C
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][3];
//			           i3=ibk(k1,2)         !k1  CA
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][2];//WHY CONNECT WITH CA?
        }




//        if(k_lineInFileIndex_inResidue==11){
//            finalFeedback.println("\nk: "+k_lineInFileIndex_inResidue+"\nk1: "+k1_residueInd+", sidelb[id1][1][id2-4]: "+
//            parameters.sideGeometry.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][1][id2_atomInResInd-4]+
//            "\nsidelb[id1][2][id2-4]: "+parameters.sideGeometry.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][2][id2_atomInResInd-4]);
//        }
//			        elseif(id2.gt.4) then
        if(id2_atomInResInd>4){//FINDS NEIGHBORING ATOMS???
//			           k1  CB = sidelb(id1,1,2)
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd]
                    [parameters.sideGeometry.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][1][id2_atomInResInd-4]];
            //THE LENGTH OF THE LINE IN SIDEGEO IS ALWAYS THE NUMBER OF ATOMS MINUS 4!!!!!!!!
            //the first 4 getAtoms have the same geometry in each residue(N,CA,C,O).

            i3_nextAtomLineA2_orPrevPrev=-1;//this will be changed in the next if statement.  Tommer 27.4.15
            if(parameters.sideGeometry.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine
                    [id1_resCodeInd][1][id2_atomInResInd-3]>0){
                //checks that we are not at the last atom in the residue(?) 27.4.15
                for(int n_atomsInd=id2_atomInResInd+1;n_atomsInd<=15;n_atomsInd++){
                    if(parameters.sideGeometry.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][1][n_atomsInd-4]==id2_atomInResInd){
                        //this checks which atom has the atom indexed by id2 as its previous atom. that is the next atom. 27.4.15
                        i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][n_atomsInd];
                        //this puts the next atom line index in place at i3.
                        n0_flagForNextFound=1;
                    }
                }
            }
            if(i3_nextAtomLineA2_orPrevPrev<0)
                i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes
                        [k1_residueInd][parameters.sideGeometry.
                        sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine
                        [id1_resCodeInd][2][id2_atomInResInd-4]];
        }
        ibb_previousAtomLine=i2_previousAtomLineA1;


//        if((k_lineInFileIndex_inResidue<13)&&(k_lineInFileIndex_inResidue>5)){
//                        finalFeedback.println("\nn0_flagNext: "+n0_flagForNextFound+", k="+k_lineInFileIndex_inResidue+
//                                "\nid1: "+id1_resCodeInd+", id2: "+id2_atomInResInd+
//                                "\ni1: "+ i1_kReplica_lineInFileIndex+", i2: "+i2_previousAtomLineA1+", i3: "+i3_nextAtomLineA2_orPrevPrev+
//                                "\nxp[i1]: "+xp_unknown[i1_kReplica_lineInFileIndex]+", yp[i1]: "+yp_unknown[i1_kReplica_lineInFileIndex]+
//                                ", zp[i1]: "+zp_unknown[i1_kReplica_lineInFileIndex]+
//                                "\nxp[i2]: "+xp_unknown[i2_previousAtomLineA1]+", yp[i2]: "+yp_unknown[i2_previousAtomLineA1]+
//                                ", zp[i2]: "+zp_unknown[i2_previousAtomLineA1]+
//                                "\nxp[i3]: "+xp_unknown[i3_nextAtomLineA2_orPrevPrev]+", yp[i3]: "+yp_unknown[i3_nextAtomLineA2_orPrevPrev]+
//                                ", zp[i3]: "+zp_unknown[i3_nextAtomLineA2_orPrevPrev]+
//                                "\n\n");
//                }


        if(n0_flagForNextFound<0){//only one bonded heavy atom
            Variables.xn1_Vz_unnormalized[1]=xp_unknown[i2_previousAtomLineA1]-xp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[2]=yp_unknown[i2_previousAtomLineA1]-yp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[3]=zp_unknown[i2_previousAtomLineA1]-zp_unknown[i1_kReplica_lineInFileIndex];
        }

        else{
            // Chen wild guess
            if (i2_previousAtomLineA1 == -1) return -1;
            Variables.xn1_Vz_unnormalized[1]=xp_unknown[i2_previousAtomLineA1]-xp_unknown[i1_kReplica_lineInFileIndex]+
                    xp_unknown[i3_nextAtomLineA2_orPrevPrev]-xp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[2]=yp_unknown[i2_previousAtomLineA1]-yp_unknown[i1_kReplica_lineInFileIndex]+
                    yp_unknown[i3_nextAtomLineA2_orPrevPrev]-yp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[3]=zp_unknown[i2_previousAtomLineA1]-zp_unknown[i1_kReplica_lineInFileIndex]+
                    zp_unknown[i3_nextAtomLineA2_orPrevPrev]-zp_unknown[i1_kReplica_lineInFileIndex];

        }

        v2_r13_unnormalized[1]=xp_unknown[i3_nextAtomLineA2_orPrevPrev]-xp_unknown[i1_kReplica_lineInFileIndex];
        v2_r13_unnormalized[2]=yp_unknown[i3_nextAtomLineA2_orPrevPrev]-yp_unknown[i1_kReplica_lineInFileIndex];
        v2_r13_unnormalized[3]=zp_unknown[i3_nextAtomLineA2_orPrevPrev]-zp_unknown[i1_kReplica_lineInFileIndex];


        v1_Vy_unnormalized=xdot_vectorMult(Variables.xn1_Vz_unnormalized,v2_r13_unnormalized);
        Variables.xd1_Vx_unnormalized=xdot_vectorMult(v1_Vy_unnormalized, Variables.xn1_Vz_unnormalized);


        double norm_xd=Math.sqrt(Variables.xd1_Vx_unnormalized[1]*Variables.xd1_Vx_unnormalized[1]+
                Variables.xd1_Vx_unnormalized[2]*Variables.xd1_Vx_unnormalized[2]+
                Variables.xd1_Vx_unnormalized[3]*Variables.xd1_Vx_unnormalized[3]);
        double norm_xn=Math.sqrt(Variables.xn1_Vz_unnormalized[1]*Variables.xn1_Vz_unnormalized[1]+
                Variables.xn1_Vz_unnormalized[2]*Variables.xn1_Vz_unnormalized[2]+
                Variables.xn1_Vz_unnormalized[3]*Variables.xn1_Vz_unnormalized[3]);

        for(int axisInd=1;axisInd<=3;axisInd++){
            Variables.xd1_Vx_unnormalized[axisInd]=Variables.xd1_Vx_unnormalized[axisInd]/norm_xd;
            Variables.xn1_Vz_unnormalized[axisInd]=Variables.xn1_Vz_unnormalized[axisInd]/norm_xn;
        }
        return ibb_previousAtomLine;

    }



    private static int caldplane_unknown_clean(double xp_unknown[], double yp_unknown[], double[] zp_unknown,
                                               int k1_residueInd, int k_lineInFileIndex_inResidue, int id1_resCodeInd, int id2_atomInResInd,
                                               int[][] ibk_idByResPlaceInFileAndAtomPlaceInRes) throws FileNotFoundException{

        int ibb_previousAtomLine; //this will be returned.


        //double[][][] sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine = new double[21][5][16];
        //this means it will be created anew each time the function is called!
        //this variable is not used in the function!

        int l_lineLengthInGeometryFile;
        double[] v1_Vy_unnormalized= new double[4], v2_r13_unnormalized=new double[4];//, v3_unknown=new double[4];

        if (Variables.init_flagFirstTimeIncaldplane){//in the first time, we need to read sidelb
            //from side_geometry.dat
            Variables.init_flagFirstTimeIncaldplane=false;
            String filename = Variables.base_goapalonePath+"/side_geometry.dat";
            String[] lineElements;
            Scanner scanner_10_sideGeometry = new Scanner (new File(filename));
            for (int i_residueNumberInGeoFile=1;i_residueNumberInGeoFile<=20;i_residueNumberInGeoFile++){
                lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                l_lineLengthInGeometryFile=Integer.parseInt(lineElements[2]);
                lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                for (int placeInGeoLine=1;placeInGeoLine<=l_lineLengthInGeometryFile;placeInGeoLine++){
                    //sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][1][placeInGeoLine]=
                    //	Double.parseDouble(lineElements[2*placeInGeoLine]);
                    Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][1][placeInGeoLine]=
                            Integer.parseInt(lineElements[2*placeInGeoLine+1]);

                }
                lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                for (int placeInGeoLine=1;placeInGeoLine<=l_lineLengthInGeometryFile;placeInGeoLine++){
                    //sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][2][placeInGeoLine]=Double.parseDouble
                    //	(lineElements[2*placeInGeoLine]);
                    Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][2][placeInGeoLine]=Integer.parseInt
                            (lineElements[2*placeInGeoLine+1]);

                }
                lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                for (int placeInGeoLine=1;placeInGeoLine<=l_lineLengthInGeometryFile;placeInGeoLine++){
                    //sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][3][placeInGeoLine]=Double.parseDouble
                    //	(lineElements[2*placeInGeoLine]);
                    Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][3][placeInGeoLine]=Integer.parseInt
                            (lineElements[2*placeInGeoLine+1]);

                }

                lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                for (int placeInGeoLine=1;placeInGeoLine<=l_lineLengthInGeometryFile;placeInGeoLine++){
                    //	sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][4][placeInGeoLine]=Double.parseDouble
                    //		(lineElements[2*placeInGeoLine]);
                    Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][3][placeInGeoLine]=Integer.parseInt
                            (lineElements[2*placeInGeoLine+1]);////THIS IS REDUNDANT!!!
                    //we have 3 again here because sidelb only has size 3 in the relevant dimension

                }//this is redundant. 4.5.15

            }
            scanner_10_sideGeometry.close();
        }
        int i1_kReplica_lineInFileIndex = k_lineInFileIndex_inResidue, i2_previousAtomLineA1=0, i3_nextAtomLineA2_orPrevPrev=0;
//			         n0=-1
        int n0_flagForNextFound=-1;
//			         if(id2.eq.1) then ! N
        if(id2_atomInResInd==1){//Nitrogen //FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1-1,3)        ! k1-1's C
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd-1][3];

//			           i3=ibk(k1,2)          ! k1's CA
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][2];
//			           n0=1
            n0_flagForNextFound=1;
        }

//			         elseif(id2.eq.2) then !  CA
        if(id2_atomInResInd==2){ //FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1,1)         !k1  N
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][1];

//			           i3=ibk(k1,3)         !k1  C
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][3];
//			           n0=1
            n0_flagForNextFound=1;
        }

//			         elseif(id2.eq.3) then !  C
        if(id2_atomInResInd==3){//FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1,2)         !k1  CA
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][2];

//			           i3=ibk(k1,4)         !k1  O
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][4];
//			           n0=1
            n0_flagForNextFound=1;
        }

//			        elseif(id2.eq.4) then !  O
        if(id2_atomInResInd==4){//FINDS NEIGHBORING ATOMS???
//			           i2=ibk(k1,3)         !k1  C
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][3];

//			           i3=ibk(k1,2)         !k1  CA
            i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][2];//WHY CONNECT WITH CA?
        }



//                if(k_lineInFileIndex_inResidue==11){
//                    finalFeedback.println("\nk: "+k_lineInFileIndex_inResidue+"\nk1: "+k1_residueInd+", sidelb[id1][1][id2-4]: "+
//                    Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][1][id2_atomInResInd-4]+
//                    "\nsidelb[id1][2][id2-4]: "+Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][2][id2_atomInResInd-4]);
//                }
//			        elseif(id2.gt.4) then
        if(id2_atomInResInd>4){//FINDS NEIGHBORING ATOMS???
//			           k1  CB = sidelb(id1,1,2)
            i2_previousAtomLineA1=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd]
                    [Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][1][id2_atomInResInd-4]];
            //THE LENGTH OF THE LINE IN SIDEGEO IS ALWAYS THE NUMBER OF ATOMS MINUS 4!!!!!!!!
            //the first 4 getAtoms have the same geometry in each residue(N,CA,C,O).

            i3_nextAtomLineA2_orPrevPrev=-1;//this will be changed in the next if statement.  Tommer 27.4.15
            if(Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine
                    [id1_resCodeInd][1][id2_atomInResInd-3]>0){
                //checks that we are not at the last atom in the residue(?) 27.4.15
                for(int n_atomsInd=id2_atomInResInd+1;n_atomsInd<=15;n_atomsInd++){
                    if(Variables.sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[id1_resCodeInd][1][n_atomsInd-4]==id2_atomInResInd){
                        //this checks which atom has the atom indexed by id2 as its previous atom. that is the next atom. 27.4.15
                        i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes[k1_residueInd][n_atomsInd];
                        //this puts the next atom line index in place at i3.
                        n0_flagForNextFound=1;
                    }
                }
            }
            if(i3_nextAtomLineA2_orPrevPrev<0)
                i3_nextAtomLineA2_orPrevPrev=ibk_idByResPlaceInFileAndAtomPlaceInRes
                        [k1_residueInd][Variables.
                        sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine
                        [id1_resCodeInd][2][id2_atomInResInd-4]];
        }
        ibb_previousAtomLine=i2_previousAtomLineA1;



//            if((k_lineInFileIndex_inResidue<13)&&(k_lineInFileIndex_inResidue>5)){
//                        finalFeedback.println("\nn0_flagNext: "+n0_flagForNextFound+", k="+k_lineInFileIndex_inResidue+
//                                "\nid1: "+id1_resCodeInd+", id2: "+id2_atomInResInd+
//                                "\ni1: "+ i1_kReplica_lineInFileIndex+", i2: "+i2_previousAtomLineA1+", i3: "+i3_nextAtomLineA2_orPrevPrev+
//                                "\nxp[i1]: "+xp_unknown[i1_kReplica_lineInFileIndex]+", yp[i1]: "+yp_unknown[i1_kReplica_lineInFileIndex]+
//                                ", zp[i1]: "+zp_unknown[i1_kReplica_lineInFileIndex]+
//                                "\nxp[i2]: "+xp_unknown[i2_previousAtomLineA1]+", yp[i2]: "+yp_unknown[i2_previousAtomLineA1]+
//                                ", zp[i2]: "+zp_unknown[i2_previousAtomLineA1]+
//                                "\nxp[i3]: "+xp_unknown[i3_nextAtomLineA2_orPrevPrev]+", yp[i3]: "+yp_unknown[i3_nextAtomLineA2_orPrevPrev]+
//                                ", zp[i3]: "+zp_unknown[i3_nextAtomLineA2_orPrevPrev]+
//                                "\n\n");
//                }

        if(n0_flagForNextFound<0){//only one bonded heavy atom
            Variables.xn1_Vz_unnormalized[1]=xp_unknown[i2_previousAtomLineA1]-xp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[2]=yp_unknown[i2_previousAtomLineA1]-yp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[3]=zp_unknown[i2_previousAtomLineA1]-zp_unknown[i1_kReplica_lineInFileIndex];
        }

        else{


            Variables.xn1_Vz_unnormalized[1]=xp_unknown[i2_previousAtomLineA1]-xp_unknown[i1_kReplica_lineInFileIndex]+
                    xp_unknown[i3_nextAtomLineA2_orPrevPrev]-xp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[2]=yp_unknown[i2_previousAtomLineA1]-yp_unknown[i1_kReplica_lineInFileIndex]+
                    yp_unknown[i3_nextAtomLineA2_orPrevPrev]-yp_unknown[i1_kReplica_lineInFileIndex];
            Variables.xn1_Vz_unnormalized[3]=zp_unknown[i2_previousAtomLineA1]-zp_unknown[i1_kReplica_lineInFileIndex]+
                    zp_unknown[i3_nextAtomLineA2_orPrevPrev]-zp_unknown[i1_kReplica_lineInFileIndex];

        }

        v2_r13_unnormalized[1]=xp_unknown[i3_nextAtomLineA2_orPrevPrev]-xp_unknown[i1_kReplica_lineInFileIndex];
        v2_r13_unnormalized[2]=yp_unknown[i3_nextAtomLineA2_orPrevPrev]-yp_unknown[i1_kReplica_lineInFileIndex];
        v2_r13_unnormalized[3]=zp_unknown[i3_nextAtomLineA2_orPrevPrev]-zp_unknown[i1_kReplica_lineInFileIndex];


        v1_Vy_unnormalized=xdot_vectorMult(Variables.xn1_Vz_unnormalized,v2_r13_unnormalized);
        Variables.xd1_Vx_unnormalized=xdot_vectorMult(v1_Vy_unnormalized, Variables.xn1_Vz_unnormalized);


        double norm_xd=Math.sqrt(Variables.xd1_Vx_unnormalized[1]*Variables.xd1_Vx_unnormalized[1]+
                Variables.xd1_Vx_unnormalized[2]*Variables.xd1_Vx_unnormalized[2]+
                Variables.xd1_Vx_unnormalized[3]*Variables.xd1_Vx_unnormalized[3]);
        double norm_xn=Math.sqrt(Variables.xn1_Vz_unnormalized[1]*Variables.xn1_Vz_unnormalized[1]+
                Variables.xn1_Vz_unnormalized[2]*Variables.xn1_Vz_unnormalized[2]+
                Variables.xn1_Vz_unnormalized[3]*Variables.xn1_Vz_unnormalized[3]);

        for(int axisInd=1;axisInd<=3;axisInd++){
            Variables.xd1_Vx_unnormalized[axisInd]=Variables.xd1_Vx_unnormalized[axisInd]/norm_xd;
            Variables.xn1_Vz_unnormalized[axisInd]=Variables.xn1_Vz_unnormalized[axisInd]/norm_xn;
        }
        return ibb_previousAtomLine;

    }


    //
    //
    //
//			c	#############################
    //subroutine xdot(x,y,z)
    private static double[] xdot_vectorMult(double[] a, double[] b){
//	         real*8 x(*),y(*),z(*)
//
        double[] aXb= new double[4];
//	         z(1)=x(2)*y(3)-x(3)*y(2)
//	         z(2)=x(3)*y(1)-x(1)*y(3)
//	         z(3)=x(1)*y(2)-x(2)*y(1)
        aXb[1]= a[2]*b[3]-a[3]*b[2];
        aXb[2]= a[3]*b[1]-a[1]*b[3];
        aXb[3]=a[1]*b[2]-a[2]*b[1];
//
//	         return
        return aXb;
//	         end
    }
    //		c	#############################

    public static void calang2_calculateAngleInDeg(double xn_vec1[], double xh1_vec2[],
                                                   double[] resultAngleAndCs){
//			       real*8 xn(3),xh1(3),xo(3),ang,cs,r1,r2,xt1,xt2,co
//			       data co/57.296/
        double  r1_sumOfSquares1=0, r2_sumOfSquares2=0, xt1_vec1component, xt2_vec2component;
        final double co_DegsInRad = 57.296;
        //
//			         cs=0.0
//			         r1=0.0
//			         r2=0.0
//			       do i=1,3
        for (int dimInd=1;dimInd<=3; dimInd++){
//			          xt1=xn(i)
//			          xt2=xh1(i)
            xt1_vec1component=xn_vec1[dimInd];
            xt2_vec2component=xh1_vec2[dimInd];
//			          cs=cs+xt1*xt2
//			          r1=r1+xt1**2
//			          r2=r2+xt2**2
            resultAngleAndCs[2]=resultAngleAndCs[2]+xt1_vec1component*xt2_vec2component;
            //caculate cs - inner product of vec1 and vec2

            r1_sumOfSquares1=r1_sumOfSquares1+xt1_vec1component*xt1_vec1component;
            r2_sumOfSquares2=r2_sumOfSquares2+xt2_vec2component*xt2_vec2component;
//			       enddo
        }
        //
//			       cs=cs/sqrt(r1*r2)
        resultAngleAndCs[2]=resultAngleAndCs[2]/Math.sqrt(r1_sumOfSquares1*r2_sumOfSquares2);
        //make cs into a cosine of the angle (dividing by norms of vec1 and vec2)
        //
//			       ang=acos(cs)*co
        resultAngleAndCs[1]=Math.acos(resultAngleAndCs[2])*co_DegsInRad;
        //we get the angle in degrees
        //
//			       return
//			       end

    }


    //			 subroutine calphi(xn,xd,xh1,ang)
    public static double calphi_signedAngleBet_xh1ANDxd_PSY(double[] xn_vec1, double[] xd_vecAdditional, double[] xh1_vec2){
//		       real*8 xn(3),xh1(3),xo(3),ang,cs,r1,r2,xt1,xt2,co
//		     &  ,xd(3),xy(3)
        double ang_phi, xo_unknown[]=new double[4], cs_cosine12=0, cs2_cosine_xo_xd, cs3_cosine_xo_xy,
                r1_sumOfSquares1=0, r2_sumOfSquares2=0,
                xt1_vec1component, xt2_vec2component, xy_vecProduct_xn_xd[]=new double[4];
//		       data co/57.296/
        final double co_DegsInRad = 57.296;

//		         cs=0.0
//		         r1=0.0
//		         r2=0.0
//		       do i=1,3
        for(int dimInd=1;dimInd<=3;dimInd++){
//		          xt1=xn(i)
//		          xt2=xh1(i)
            xt1_vec1component=xn_vec1[dimInd];
            xt2_vec2component=xh1_vec2[dimInd];
//		          cs=cs+xt1*xt2
//		          r1=r1+xt1**2
//		          r2=r2+xt2**2
            cs_cosine12=cs_cosine12+xt1_vec1component*xt2_vec2component;
            r1_sumOfSquares1=r1_sumOfSquares1+xt1_vec1component*xt1_vec1component;
            r2_sumOfSquares2=r2_sumOfSquares2+xt2_vec2component*xt2_vec2component;
//		       enddo
        }

//		       cs=cs/sqrt(r1*r2)
        cs_cosine12=cs_cosine12/Math.sqrt(r1_sumOfSquares1*r2_sumOfSquares2);
        // make cs into cosine (dividing by norms)

//		       ss=sqrt(1-cs*cs+0.0001)//redundant
        //double ss_sineAbsVal=Math.sqrt(1-cs_cosine12*cs_cosine12+0.0001);//redundant...?
//		        

//				xo(1)=xh1(1)-cs*xn(1)
//		        xo(2)=xh1(2)-cs*xn(2)
//		        xo(3)=xh1(3)-cs*xn(3)
        xo_unknown[1]=xh1_vec2[1]-cs_cosine12*xn_vec1[1];
        xo_unknown[2]=xh1_vec2[2]-cs_cosine12*xn_vec1[2];
        xo_unknown[3]=xh1_vec2[3]-cs_cosine12*xn_vec1[3];

        //this is very unclear. i can't grasp the geometrical significance of this
        //vector (XO). it furthermore seems to be the case that xn and xd are always
        //orthogonal, and therefore xo cdot xd is the same as xh1 (or vt) cdot xd.
        //this needs CHECKING! 1.6.15
//
//		       cs2=xo(1)*xd(1)+xo(2)*xd(2)+xo(3)*xd(3)
//		       cs2=cs2/sqrt(xo(1)**2+xo(2)**2+xo(3)**2)
        cs2_cosine_xo_xd=xo_unknown[1]*xd_vecAdditional[1]+xo_unknown[2]*xd_vecAdditional[2]+
                xo_unknown[3]*xd_vecAdditional[3];
        cs2_cosine_xo_xd=cs2_cosine_xo_xd/Math.sqrt(xo_unknown[1]*xo_unknown[1]+
                xo_unknown[2]*xo_unknown[2]+xo_unknown[3]*xo_unknown[3]);

//		       call xdot(xn,xd,xy)
        xy_vecProduct_xn_xd=xdot_vectorMult(xn_vec1, xd_vecAdditional);

//		       cs3=xy(1)*xo(1)+xy(2)*xo(2)+xy(3)*xo(3)
        cs3_cosine_xo_xy=xy_vecProduct_xn_xd[1]*xo_unknown[1]+xy_vecProduct_xn_xd[2]*xo_unknown[2]+
                xy_vecProduct_xn_xd[3]*xo_unknown[3];

//		       if(cs2.lt.-1.0) cs2=-0.9999
        if(cs2_cosine_xo_xd<-1.0) cs2_cosine_xo_xd=-0.9999;
//		       if(cs2.gt.1.0) cs2=0.9999
        if(cs2_cosine_xo_xd>1.0) cs2_cosine_xo_xd=0.9999;

//		       ang=acos(cs2)*co
        ang_phi=Math.acos(cs2_cosine_xo_xd)*co_DegsInRad;//gives PHI in degrees. but what angle is phi??
//		       if(cs3.lt.0) ang=-ang
        if(cs3_cosine_xo_xy<0) ang_phi=-ang_phi;
//
//		       return


        //-----------------------------------------checkup 1.6.15-----------------------------
				/*outwriter.println("\ncheckup in calphi 1.6.15");
				outwriter.println("xn_vec1: "+xn_vec1[1]+", "+xn_vec1[2]+", "+xn_vec1[3]);
				outwriter.println("xh1_vec2: "+xh1_vec2[1]+", "+xh1_vec2[2]+", "+xh1_vec2[3]);
				outwriter.println("cs12: "+ cs_cosine12);
				outwriter.println("xo: "+xo_unknown[1]+", "+xo_unknown[2]+", "+xo_unknown[3]);
				outwriter.println("cs2: "+cs2_cosine_xo_xd);*/

        return ang_phi;
//		       end
    }


    public static double dihedral_calc_Xhi(double[] x_unknown, double[] y_unknown, double[] z_unknown){
//  			subroutine dihedral (x,y,z,dih)
//		      parameter (Nn=4)
        final int Nn_dimension=4;
        final double PI=Math.PI;
        int ipi_unknown=1, jpi_unknown=2,kpi_unknown=3,lpi_unknown=4;
//		      real*8  x(Nn),y(Nn), z(Nn),dih
//		      real*8 mconst
        double dih_unknown, mconst_unknown, dihedral_angle=0;
//		c
//		c bond angle distribution
//		c
//		      pi = 4.d0*datan(1.d0)
//		         i  = 1
//		         ipi = i
//		         jpi = i+1
//		         kpi = i+2
//		         lpi = i+3
//		c
//		c dihedral angle
//		         dihedral_angle = 0.d0
//		c
//		         rijx = x(ipi) - x(jpi)
//		         rijy = y(ipi) - y(jpi)
//		         rijz = z(ipi) - z(jpi)
//		         rjkx = x(jpi) - x(kpi)
//		         rjky = y(jpi) - y(kpi)
//		         rjkz = z(jpi) - z(kpi)
//		         rklx = x(kpi) - x(lpi)
//		         rkly = y(kpi) - y(lpi)
//		         rklz = z(kpi) - z(lpi)

        double rijx_,rijy_,rijz_,rjkx_,rjky_,rjkz_,rklx_,rkly_,rklz_,
                ax_, ay_, az_, bx_, by_, bz_;
        rijx_=x_unknown[ipi_unknown]-x_unknown[jpi_unknown];
        rijy_=y_unknown[ipi_unknown]-y_unknown[jpi_unknown];
        rijz_=z_unknown[ipi_unknown]-z_unknown[jpi_unknown];
        //vector from ipi to jpi

        rjkx_=x_unknown[jpi_unknown]-x_unknown[kpi_unknown];
        rjky_=y_unknown[jpi_unknown]-y_unknown[kpi_unknown];
        rjkz_=z_unknown[jpi_unknown]-z_unknown[kpi_unknown];
        //vector from jpi to kpi

        rklx_=x_unknown[kpi_unknown]-x_unknown[lpi_unknown];
        rkly_=y_unknown[kpi_unknown]-y_unknown[lpi_unknown];
        rklz_=z_unknown[kpi_unknown]-z_unknown[lpi_unknown];
        //vector from kpi to lpi


//		c
//		         ax = rijy*rjkz - rijz*rjky
//		         ay = rijz*rjkx - rijx*rjkz
//		         az = rijx*rjky - rijy*rjkx
//		         bx = rjky*rklz - rjkz*rkly
//		         by = rjkz*rklx - rjkx*rklz
//		         bz = rjkx*rkly - rjky*rklx
        ax_= rijy_*rjkz_-rijz_*rjky_;
        ay_=rijz_*rjkx_-rijx_*rjkz_;
        az_=rijx_*rjky_-rijy_*rjkx_;
        //a is the vector product betwween the two first distance vectors

        bx_= rjky_*rklz_-rjkz_*rkly_;
        by_=rjkz_*rklx_-rjkx_*rklz_;
        bz_=rjkx_*rkly_-rjky_*rklx_;
        //b is the vector product betwween the two last distance vectors

        //outWriter.println("\ncheckup 3.6.15 \na: "+ax_+", "+ay_+", "+az_);//checkup
        //outWriter.println("b: "+bx_+", "+by_+", "+bz_);//checkup

//		C
//		C     set to MCONST if smaller than MCONST
//		         mconst=1.0d-10
        mconst_unknown=1.0/Math.pow(10, 10);
//		         RG=SQRT(MAX(MCONST,RJKX*RJKX+RJKY*RJKY+RJKZ*RJKZ))
        double rg_RjkNorm_withLowerLimit=Math.sqrt(Math.max(mconst_unknown, rjkx_*rjkx_+rjky_*rjky_+rjkz_*rjkz_));
        //rg is set to the norm of Rjk or to 10^-5, whichever is bigger.

//		         RGR=1.d0/RG
        //double rgr_RjkNormReciprocal_withUpperLimit=1/rg_RjkNorm_withLowerLimit;
//		         RA2R=1.d0/MAX(MCONST,AX*AX+AY*AY+AZ*AZ)
        double ra2r_aNormSqrReciprocal_withUpperLimit=1.0/Math.max(mconst_unknown, ax_*ax_+ay_*ay_+az_*az_);
//		         RB2R=1.d0/MAX(MCONST,BX*BX+BY*BY+BZ*BZ)
        double rb2r_bNormSqrReciprocal_withUpperLimit=1.0/Math.max(mconst_unknown, bx_*bx_+by_*by_+bz_*bz_);
//		         RABR=SQRT(RA2R*RB2R)
        double rabr_normProductABreciprocal=Math.sqrt(ra2r_aNormSqrReciprocal_withUpperLimit*
                rb2r_bNormSqrReciprocal_withUpperLimit);
        //outWriter.println("rg: "+rg_RjkNorm_withLowerLimit+", rgr: "+rgr_RjkNormReciprocal_withUpperLimit+
        //	", ra2r: "+ra2r_aNormSqrReciprocal_withUpperLimit+", rb2r: "+rb2r_bNormSqrReciprocal_withUpperLimit+
        //", rabr: "+rabr_normProductABreciprocal);//checkup
//		CP=COS(PHI)
//		         CP=RABR*(AX*BX+AY*BY+AZ*BZ)
        double cp_cosineAB=rabr_normProductABreciprocal*(ax_*bx_+ay_*by_+az_*bz_);
//		C SP=SIN(PHI)
//		C which can be simplify to sin(phi)=|G|H.A/(|A|.|B|)
//		Cab...B950603 06/29/95
//		Cab        SP=RG*RABR*(AX*HX+AY*HY+AZ*HZ)
//		         SP=-RG*RABR*(AX*RKLX+AY*RKLY+AZ*RKLZ)
        double sp_parallelopipedVolumeDividedByNorms=-rg_RjkNorm_withLowerLimit*rabr_normProductABreciprocal*
                (ax_*rklx_+ay_*rkly_+az_*rklz_);
//		Cab...
//		C
//		         if(cp.gt.1.d0) then
        double theta_angleBetweenNormsToPlanes;
				/*if (printCounter==11172){//checkup 30.6.15
					System.out.println("checkup 30.6\ncp_cosineAB: "+cp_cosineAB+"\n\n");
				}*/
        if(cp_cosineAB>1){
//		c            write(99,*) 'dh:', cp, ipi,jpi,kpi,lpi
//		            cp = 1.d0
//		            theta  = 0.d0
            theta_angleBetweenNormsToPlanes=0;
            cp_cosineAB=1;
        }
        //		         else if(cp.lt.-1.d0) then
        else if(cp_cosineAB<-1){
            //		c            write(99,*) 'dh:', cp, ipi,jpi,kpi,lpi
            //		            cp = -1.d0
            cp_cosineAB=-1;
            //		            theta  = 180.d0
            theta_angleBetweenNormsToPlanes=180;
        }
        else{
            //		         else
            //		            theta    = acos(cp)
            theta_angleBetweenNormsToPlanes=Math.acos(cp_cosineAB);
            //		            theta    = theta*180/pi
            theta_angleBetweenNormsToPlanes=theta_angleBetweenNormsToPlanes*180/PI;
            //		         endif
        }
//		         if(sp.lt.0.0) theta = -theta          
//		         dih = theta
        if (sp_parallelopipedVolumeDividedByNorms<0) theta_angleBetweenNormsToPlanes=-theta_angleBetweenNormsToPlanes;
        //dih_unknown=theta_angleBetweenNormsToPlanes;
//		         return
        return theta_angleBetweenNormsToPlanes;
        //return dih_unknown;
//		      end
    }


}
