package meshi.energy.goap;

import meshi.energy.AbstractEnergy;
import meshi.energy.EnergyCreator;
import meshi.geometry.DistanceMatrix;
import meshi.molecularElements.Protein;
import meshi.util.CommandList;
import meshi.util.KeyWords;
import meshi.util.Utils;
import meshi.util.file.MeshiLineReader;
import meshi.util.info.InfoType;

import java.io.*;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

/**
 * Created by chen on 05/07/2015.
 */
public class GoapCreator extends EnergyCreator{
    public static final int ibin_unknown=20; //this is originaly declared inside a "pararmeter(...)"
    int[] map_unknown = new int[501];
    public static  int ibinme_unknown;
    public static  int mapnum_numOfLinesToBeRead;
    String lineFromFortReader;
    final double dh_unknown = 2./ibin_unknown;
    public  String parametersPath;
    String[] lineElements = new String[22];
    Fort21 fort21;
    Fort31 fort31;
    Charges charges;
    SideGeometry sideGeometry;
    private GoapParameters parameters = null;
    public GoapCreator() {
        super(InfoType.GOAP);
    }

    public AbstractEnergy createEnergyTerm(Protein protein,DistanceMatrix distanceMatrix, CommandList commands) {
        if (parameters ==  null) {
            parametersPath = commands.firstWord(KeyWords.PARAMETERS_DIRECTORY).secondWord() + "/meshiPotential/GOAP/";
            Utils.println("checkup 13.7.15, parametersPath: " + parametersPath);
            try {
                //fort21 = getFort21();// siditom change 15.2.2018 - read object from memory instead of text file.
                Utils.println("checkup 15.2.18, fort21");
                fort21 = getFort21Object();
            } catch (Exception ex) {
                throw new RuntimeException(ex);
            }

            try {
                Utils.println("checkup 15.2.18, fort31");
                //fort31 = getFort31(); // siditom change 15.2.2018 - read object from memory instead of text file.
                fort31 = getFort31Object();
            } catch (Exception ex) {
                throw new RuntimeException(ex);
            }

            try {
                charges = getCharges();
            } catch (IOException ex) {
                throw new RuntimeException(ex);
            }

            // try {
            sideGeometry = getSideGeometry();
            //}   catch ( IOException ex) {throw new RuntimeException(ex);}

            parameters = new GoapParameters(fort21, fort31, charges, sideGeometry);
        }
        return term = new Goap(protein,parameters, new GoapInfo());
    }

    public  Charges getCharges() throws IOException{
        //
//			open(unit=10,file=base(1:ips)//'/charge_inp.dat',
//		     &                     status='old')

        MeshiLineReader data_from_chargeReader_10=
                new MeshiLineReader(parametersPath+"charge_inp.dat");
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

        return new Charges(cind_atomsByOrder,chg_chargesByOrder, resn_threeLetterResidueCodes,ianum_numbersOfAtomsInResidues);
    }

    public  Fort31 getFort31Object() throws Exception {
        FileInputStream fis = new FileInputStream(parametersPath+"fort31_db");
        GZIPInputStream gs = new GZIPInputStream(fis);
        ObjectInputStream ois = new ObjectInputStream(gs);
        Fort31 fort31 = (Fort31) ois.readObject();
        ois.close();
        fis.close();
        return fort31;
    }

    public  Fort31 getFort31() throws IOException {
        FortranArray7Dim_cnttheta cnttheta_unknown = new FortranArray7Dim_cnttheta();
        //		        open(unit=21,file=base(1:ips)//'/fort.31_g72_noshift5_new',
//		     &               status='old')
        MeshiLineReader data_from_fort31Reader_21 =
                new MeshiLineReader(parametersPath+"fort.31_g72_noshift5_new");
//		         read(21,*)
        data_from_fort31Reader_21.readLine();//reads comment line
//		         read(21,*) ibinme,mapnum,ig
        Scanner scanner = new Scanner(data_from_fort31Reader_21.readLine());
        ibinme_unknown = scanner.nextInt();
        mapnum_numOfLinesToBeRead = scanner.nextInt();
        int ig_s_parameter = scanner.nextInt();
        scanner.close();
//		         if(ibinme.gt.ibin) then
//		         write(*,*) 'bin exceeds limit:',ibin,'!'
//		         stop
//		         endif
        if (ibinme_unknown > ibin_unknown) throw new RuntimeException("bin exceeds limit: " + ibin_unknown + "!");
//		         do i=1,mapnum
//		         read(21,*) map(i)
//		c tommer this runs over all the previously read values in map() from the previous file. they might very well have the exact same value (should check!)
//		         enddo
        for (int i = 1; i <= mapnum_numOfLinesToBeRead; i++) {//this collects again values for map[].  (same values) 13.3.15
            scanner = new Scanner(data_from_fort31Reader_21.readLine());
            map_unknown[i] = scanner.nextInt();
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

        lineFromFortReader = data_from_fort31Reader_21.readLine();
        float[] pairPotentialX_ppx = new float[13];
        float pairPotentialY_ppy;
        int numberOfLineBeingRead = 0;
        //long startTime = System.currentTimeMillis();
        int m, i, j, k, l, mx;
        int internalLineNumber, lx, placeInLine;
        float pairPotential_pp;
        data_from_fort31Reader_21.close();//for using File
        File fileFort31 = new File(parametersPath+ "/fort.31_g72_noshift5_new");
        scanner = new Scanner(new FileReader(fileFort31));
        for (int lineToSkipNum = 1; lineToSkipNum <= 52; lineToSkipNum++) scanner.nextLine();
        String line = scanner.nextLine();
        //String[] lineElements = new String[22];


        while (line != null) {
            boolean flag_noLine = false;
            numberOfLineBeingRead++;

            lineElements = line.split("\\s+");

            try {
                m = Integer.parseInt(lineElements[4]);
            } catch (Exception e) {
                System.out.println("line count: " + numberOfLineBeingRead);
                for (int elementNum = 0; elementNum < lineElements.length; elementNum++)
                    System.out.println("element number " + elementNum + " : " + lineElements[elementNum]);
                throw new  RuntimeException(e);

            }
            pairPotential_pp = Float.parseFloat(lineElements[5]);
            try {
                i = Integer.parseInt(lineElements[6]);
            } catch (Exception e) {
                System.out.println("line count: " + numberOfLineBeingRead);
                for (int elementNum = 0; elementNum < lineElements.length; elementNum++)
                    System.out.println("element number " + elementNum + " : " + lineElements[elementNum]);
                throw new RuntimeException(e);

            }
            //i=Integer.parseInt(lineElements[6]);
            j = Integer.parseInt(lineElements[7]);
            k = Integer.parseInt(lineElements[8]);
            l = Integer.parseInt(lineElements[9]);

            if (scanner.hasNextLine()) line = scanner.nextLine();
            else {
               Utils.println("line: " + numberOfLineBeingRead);
                break;
            }
            for (internalLineNumber = 1; (internalLineNumber <= 5) && (line != null); internalLineNumber++) {
                numberOfLineBeingRead++;
                lineElements = line.split("\\s+");
                mx = Integer.parseInt(lineElements[5]);
                for (placeInLine = 1; placeInLine <= 12; placeInLine++) {
                    pairPotentialX_ppx[placeInLine] = Float.parseFloat(lineElements[5 + placeInLine]);
                }
                pairPotentialY_ppy = Float.parseFloat(lineElements[18]);
                //INDICES have been changed (23.3.15). must be CHECKED!


                if (!((i == 21) && (k == 21))) {//is this necessary? tommer 13.3.15
                    for (lx = 1; lx <= 12; lx++) {
                        cnttheta_unknown.set_cnttheta(i, j, k, l, m, lx, mx, pairPotentialX_ppx[lx]);
                    }
                }
                if (scanner.hasNextLine()) line = scanner.nextLine();
                else {
                    Utils.println("line: " + numberOfLineBeingRead + " , inside loop. inside index: " + internalLineNumber);
                    flag_noLine = true;
                    break;
                }
            }
            if (flag_noLine) {
                if (internalLineNumber != 5) throw new RuntimeException("file ends with wrong line pattern (not 1+5)");
                break;
            }
            for (lx = 1; lx <= 12; lx++) {//symmetries
                cnttheta_unknown.set_cnttheta(k, l, i, j, m, lx, 1, cnttheta_unknown.get_cnttheta(i, j, k, l, m, lx, 3));
                cnttheta_unknown.set_cnttheta(k, l, i, j, m, lx, 2, cnttheta_unknown.get_cnttheta(i, j, k, l, m, lx, 4));
                cnttheta_unknown.set_cnttheta(k, l, i, j, m, lx, 3, cnttheta_unknown.get_cnttheta(i, j, k, l, m, lx, 1));
                cnttheta_unknown.set_cnttheta(k, l, i, j, m, lx, 4, cnttheta_unknown.get_cnttheta(i, j, k, l, m, lx, 2));
                cnttheta_unknown.set_cnttheta(k, l, i, j, m, lx, 5, cnttheta_unknown.get_cnttheta(i, j, k, l, m, lx, 5));

            }


        }





        return new Fort31(cnttheta_unknown,map_unknown,ig_s_parameter);
    }

    public Fort21 getFort21Object() throws Exception{
        FileInputStream fis = new FileInputStream(parametersPath+"fort21_db");
        GZIPInputStream gs = new GZIPInputStream(fis);
        ObjectInputStream ois = new ObjectInputStream(gs);
        Fort21 fort21 = (Fort21) ois.readObject();
        ois.close();
        fis.close();
        return fort21;
    }

    public Fort21 getFort21() throws IOException{



        MeshiLineReader data_from_fort21Reader_20 = new MeshiLineReader(parametersPath+"fort.21_1.61_2");
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





        FortranArray5Dim_pot pot_potential = new FortranArray5Dim_pot();

        lineFromFortReader=data_from_fort21Reader_20.readLine();
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
        return new Fort21(map_unknown, pot_potential);
    }



    private SideGeometry getSideGeometry()/*throws FileNotFoundException*/ {
        int l_lineLengthInGeometryFile;
        int[][][] sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine = new int[21][4][16];
            //from side_geometry.dat
            String filename = parametersPath+"side_geometry.dat";
           Utils.println("checkup 13.7.15, filename: "+filename);
            String[] lineElements;


            try {
                Scanner scanner_10_sideGeometry = new Scanner(new File(filename));
                for (int i_residueNumberInGeoFile = 1; i_residueNumberInGeoFile <= 20; i_residueNumberInGeoFile++) {
                    lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                    l_lineLengthInGeometryFile = Integer.parseInt(lineElements[2]);
                    lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                    for (int placeInGeoLine = 1; placeInGeoLine <= l_lineLengthInGeometryFile; placeInGeoLine++) {
                        //sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][1][placeInGeoLine]=
                        //	Double.parseDouble(lineElements[2*placeInGeoLine]);
                        sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][1][placeInGeoLine] =
                                Integer.parseInt(lineElements[2 * placeInGeoLine + 1]);

                    }
                    lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                    for (int placeInGeoLine = 1; placeInGeoLine <= l_lineLengthInGeometryFile; placeInGeoLine++) {
                        //sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][2][placeInGeoLine]=Double.parseDouble
                        //	(lineElements[2*placeInGeoLine]);
                        sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][2][placeInGeoLine] = Integer.parseInt
                                (lineElements[2 * placeInGeoLine + 1]);

                    }
                    lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                    for (int placeInGeoLine = 1; placeInGeoLine <= l_lineLengthInGeometryFile; placeInGeoLine++) {
                        //sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][3][placeInGeoLine]=Double.parseDouble
                        //	(lineElements[2*placeInGeoLine]);
                        sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][3][placeInGeoLine] = Integer.parseInt
                                (lineElements[2 * placeInGeoLine + 1]);

                    }

                    lineElements = scanner_10_sideGeometry.nextLine().split("\\s+");
                    for (int placeInGeoLine = 1; placeInGeoLine <= l_lineLengthInGeometryFile; placeInGeoLine++) {
                        //	sidegeo_doubleSideData_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][4][placeInGeoLine]=Double.parseDouble
                        //		(lineElements[2*placeInGeoLine]);
                        sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine[i_residueNumberInGeoFile][3][placeInGeoLine] = Integer.parseInt
                                (lineElements[2 * placeInGeoLine + 1]);////THIS IS REDUNDANT!!!
                        //we have 3 again here because sidelb only has size 3 in the relevant dimension

                    }//this is redundant. 4.5.15

                }
                scanner_10_sideGeometry.close();
                return new SideGeometry(sidelb_intSideDataByResinFileAndLineinResAndPlaceInLine);

            }
            catch(IOException ex){
                ex.printStackTrace();
                throw new RuntimeException("can't find sideGeometry file. Checkup Tommer 13.7.15");
            }
    }


}
