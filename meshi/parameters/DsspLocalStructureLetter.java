package meshi.parameters;


/**
 * Created by user on 14/02/2017.
 */
public enum DsspLocalStructureLetter {

    HELIX("H"),         //0
    SHEET("E"),         //1
    BETABRIDGE("B"),    //2
    THREE10HELIX("G"),  //3
    PIHELIX("I"),       //4
    BEND("S"),          //5
    TURN("T"),          //6
    HB_OandN("X"),      //7
    HB_onlyO(">"),      //8
    HB_onlyN("<"),      //9
    BRACKET_3RES("3"),  //10
    BRACKET_4RES("4"),  //11
    BRACKET_5RES("5"),  //12
    CHIRALITY_NEGATIVE("-"), //13
    CHIRALITY_POSITIVE("+"), //14
    PARALLEL("a"),       //15
    ANTI_PARALLEL("A"),  //16
    GAP("_"),            //17
    UNK("!"),            //18

    // Reduction letters - additional semantic letters
    //STR2 additional letters
    TWO_SIDED_PARALLEL_BETA("P"), //19
    TWO_SIDED_ANTI_PARALLEL_BETA("A"), //20
    TWO_SIDED_MIXED_BETA("M"), //21
    ONE_SIDED_PARALLEL_BETA("Q"), //22
    ONE_SIDED_ANTI_PARALLEL_BETA("Z"), //23

    //DSSP3 additional letters
    COIL("C"); //24

    public String getNameOneLetter() {
        return nameOneLetter;
    }

    private final String nameOneLetter;
    private DsspLocalStructureLetter(String nameOneLetter) {
        this.nameOneLetter = nameOneLetter;
    }



    public static DsspLocalStructureLetter dsspLocalStructureLetter(int loc, char c) {
        char newC = Character.toUpperCase(c);
        DsspLocalStructureLetter dls = null;
        if (loc == 7 || loc == 8) {
            if (Character.isLowerCase(c)) dls = DsspLocalStructureLetter.PARALLEL;
            else if (Character.isUpperCase(c)) dls = DsspLocalStructureLetter.ANTI_PARALLEL;
            else dls = DsspLocalStructureLetter.GAP;
        } else { //if
            for (DsspLocalStructureLetter sec_letter : DsspLocalStructureLetter.values()) {
                if (sec_letter != DsspLocalStructureLetter.PARALLEL &&
                        sec_letter != DsspLocalStructureLetter.ANTI_PARALLEL &&
                        sec_letter.nameOneLetter.charAt(0) == newC)
                    dls = sec_letter;
            }//for
        }//else
        if ( dls!=null && dls.isLegalbyLocation(loc))
            return dls;
        else throw new RuntimeException("Undefined Dssp Local Structure Letter for location " +loc+": " + newC + "\n" + "Please take a look at meshi.parameters.DsspLocalStructureLetter.");
    }

    //return true if the letter can be on the first position of the dssp structure definition, else - false
    public boolean isLegalbyLocation(int loc){
        if (    (loc == 1 && isLoc1()) ||
                (loc == 2 && isLoc2_4()) || (loc == 3 && isLoc2_4()) || (loc == 4 && isLoc2_4()) ||
                (loc == 5 && isLoc5()) ||
                (loc ==6 && isLoc6()) ||
                (loc ==7 && isLoc7_8()) ||
                (loc ==8 && isLoc7_8())
                )
            return true;
        else return false;
    }
    public boolean isLoc1(){
        if (this == HELIX ||
                this == SHEET ||
                this == BETABRIDGE ||
                this == THREE10HELIX ||
                this == PIHELIX ||
                this == BEND ||
                this == TURN ||
                this == GAP)
            return true;
        else return false;
    }

    public boolean isLoc2_4(){
        if (this == HB_OandN ||
                this == HB_onlyO ||
                this == HB_onlyN ||
                this == BRACKET_3RES ||
                this == BRACKET_4RES ||
                this == BRACKET_5RES ||
                this == GAP) return true;
        else return false;
    }
    public boolean isLoc5(){
        if (this == BEND || this == GAP)
            return true;
        else return false;
    }
    public boolean isLoc6(){
        if (this == CHIRALITY_NEGATIVE ||
                this == CHIRALITY_POSITIVE ||
                this == GAP) return true;
        else return false;
    }
    public boolean isLoc7_8(){
        if (this == PARALLEL ||
                this == ANTI_PARALLEL ||
                this == GAP) return true;
        else return false;
    }


}
