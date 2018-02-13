package meshi.parameters;

import meshi.util.MeshiAttribute;

/**
 * Created by siditom on 14/02/2017.
 *
 * Number of possible Dssp structure full strings is: 49,392
 */
public class LocalStructure implements MeshiAttribute, Comparable<LocalStructure> {



    private DsspLocalStructureLetter[] dls;
    private final int structureLength = 8;
    private boolean isAssigned = false;

    //UnAssigned LocalStructure instance.
    public LocalStructure(){
        dls = new DsspLocalStructureLetter[structureLength];
        dls[0]=DsspLocalStructureLetter.UNK;
        for (int i=1; i < structureLength; i++){
            dls[i] = DsspLocalStructureLetter.dsspLocalStructureLetter(i + 1,'_');
        }
        isAssigned = false;
    }

    public LocalStructure(LocalStructure ls){
        dls = new DsspLocalStructureLetter[structureLength];
        for (int i=0; i < structureLength; i++){
            dls[i] = ls.getLocalStructureLetter(i);
        }
        isAssigned = true;
    }

    public LocalStructure(String dsspResidueStructureInfo){
        this(dsspResidueStructureInfo.toCharArray());
    }

    public LocalStructure(char[] dsspResidueStructureInfo){
        if (dsspResidueStructureInfo == null || dsspResidueStructureInfo.length != structureLength) throw new RuntimeException("Dssp Residue Structure Info is not enough: "+ dsspResidueStructureInfo.length +"\n" + "Please take a look at meshi.parameters.LocalStructure.");
        dls = new DsspLocalStructureLetter[8];
        for (int i=0; i < structureLength; i++){
            dls[i] = DsspLocalStructureLetter.dsspLocalStructureLetter(i + 1,dsspResidueStructureInfo[i]);
        }
        isAssigned = true;
    }



    public DsspLocalStructureLetter getLocalStructureLetter(int iLetter) {
        return dls[iLetter];
    }

    @Override
    public String toString(){
        String str = "";
        for (int i=0; i < structureLength; i++){
            str += dls[i].getNameOneLetter();
        }
        return str;
    }


    public boolean equals(LocalStructure other, LocalStructureAlphabetType type){
        boolean isEqual = false;
        switch (type){
            case DSSP3:
                isEqual = (this.getDSSP3Reduction() == other.getDSSP3Reduction());
                break;
            case DSSP7:
                isEqual = (this.getDSSP7Reduction() == other.getDSSP7Reduction());
                break;
            case STR2:
                isEqual = (this.getSTR2Reduction() == other.getSTR2Reduction());
                break;
            case DSSP33:
                isEqual = (this.compareTo(other)==0);
                if (!isEqual && !LocalStructure.isDSSP30Structure(this) && !LocalStructure.isDSSP30Structure(other)) isEqual = (this.getDSSP3Reduction() == this.getDSSP3Reduction());
                break;
            case DSSP103:
                isEqual = (this.compareTo(other)==0);
                if (!isEqual && !LocalStructure.isDSSP100Structure(this) && !LocalStructure.isDSSP100Structure(other)) isEqual = (this.getDSSP3Reduction() == this.getDSSP3Reduction());
                break;
                //System.out.println("***** LocalStructureAlphabetType - NOT SUPPPORTED - for more information visit: meshi.parameters.LocalStructure *****");

            default:
                throw new RuntimeException("Unsupported LocalStructureAlphabetType: " + type + "\n" + "Please take a look at meshi.parameters.LocalStructure.");
        }
        return isEqual;
    }
    public DsspLocalStructureLetter getDSSPReduction(LocalStructureAlphabetType localStructureAlphabetType){
        if (localStructureAlphabetType == LocalStructureAlphabetType.DSSP3) return this.getDSSP3Reduction();
        else if (localStructureAlphabetType == LocalStructureAlphabetType.DSSP7) return this.getDSSP7Reduction();
        else if (localStructureAlphabetType == LocalStructureAlphabetType.STR2) return this.getSTR2Reduction();
        else throw new RuntimeException("LocalStructure - getDSSPReduction - Error - Unknown reduction for LocalStructureAlphabetType "+ localStructureAlphabetType+ ".");

    }
    public DsspLocalStructureLetter getDSSP7Reduction(){
        if (this.dls[0] == DsspLocalStructureLetter.GAP) {
            return DsspLocalStructureLetter.COIL;
        }
        return this.dls[0];
    }
    public DsspLocalStructureLetter getDSSP3Reduction(){
        DsspLocalStructureLetter firstLetter = this.dls[0];
        DsspLocalStructureLetter newLetter;
        switch (firstLetter){
            case THREE10HELIX:
            case PIHELIX:
            case HELIX:
                newLetter = DsspLocalStructureLetter.HELIX;
                break;
            case BETABRIDGE:
            case SHEET:
                newLetter = DsspLocalStructureLetter.SHEET;
                break;
            case UNK:
                newLetter = DsspLocalStructureLetter.UNK;
                break;
            case BEND:
            case TURN:
            case GAP:
            default:
                newLetter = DsspLocalStructureLetter.COIL;
                break;
        }
        return newLetter;
    }

    public DsspLocalStructureLetter getSTR2Reduction(){
        DsspLocalStructureLetter firstLetter = this.dls[0];
        if (firstLetter != DsspLocalStructureLetter.SHEET)
            return this.getDSSP7Reduction();

        DsspLocalStructureLetter bp1 = this.dls[6];
        DsspLocalStructureLetter bp2 = this.dls[7];

        if ( (bp1 == DsspLocalStructureLetter.GAP) && (bp2 == DsspLocalStructureLetter.GAP)) return DsspLocalStructureLetter.SHEET;
        else if ( (bp1 == DsspLocalStructureLetter.PARALLEL) && (bp2 == DsspLocalStructureLetter.PARALLEL)) return DsspLocalStructureLetter.TWO_SIDED_PARALLEL_BETA;
        else if ( (bp1 == DsspLocalStructureLetter.PARALLEL) && (bp2 == DsspLocalStructureLetter.ANTI_PARALLEL)) return DsspLocalStructureLetter.TWO_SIDED_MIXED_BETA;
        else if ( (bp1 == DsspLocalStructureLetter.ANTI_PARALLEL) && (bp2 == DsspLocalStructureLetter.ANTI_PARALLEL)) return DsspLocalStructureLetter.TWO_SIDED_ANTI_PARALLEL_BETA;
        else if ( (bp1 == DsspLocalStructureLetter.ANTI_PARALLEL) && (bp2 == DsspLocalStructureLetter.GAP)) return DsspLocalStructureLetter.ONE_SIDED_ANTI_PARALLEL_BETA;
        else if ( (bp1 == DsspLocalStructureLetter.GAP) && (bp2 == DsspLocalStructureLetter.ANTI_PARALLEL)) return DsspLocalStructureLetter.ONE_SIDED_ANTI_PARALLEL_BETA;
        else if ( (bp1 == DsspLocalStructureLetter.PARALLEL) && (bp2 == DsspLocalStructureLetter.GAP)) return DsspLocalStructureLetter.ONE_SIDED_PARALLEL_BETA;
        else if ( (bp1 == DsspLocalStructureLetter.GAP) && (bp2 == DsspLocalStructureLetter.PARALLEL)) return DsspLocalStructureLetter.ONE_SIDED_PARALLEL_BETA;

        return DsspLocalStructureLetter.SHEET;
    }


    public String getStructureFullString(){
        String ans = "";
        for (int i=0;i<structureLength; i++){
            ans += dls[i].getNameOneLetter();
        }
        return ans;
    }

    @Override
    public int compareTo(LocalStructure o1) {
        return this.getStructureFullString().compareTo(o1.getStructureFullString());
    }

    private static boolean isDSSP30Structure(LocalStructure ls){
        return     LocalStructure.isDSSPxStructure(ls,30);
    }
    private static boolean isDSSP100Structure(LocalStructure ls){
        return     LocalStructure.isDSSPxStructure(ls,100);
    }

    private static boolean isDSSPxStructure(LocalStructure ls,int x){
        for (int iStructure=0; iStructure<30; iStructure++){
            if (ls.compareTo(new LocalStructure(LocalStructure.DsspSTRUCTURESortedByFrequency[iStructure]))==0) return true;
        }
        return false;
    }



    public int key() {
        return SS_DSSP;
    }
    //private final static LocalStructure[] Dssp30MostFrequentStructures = {};

    //DSSP structure columns possibilities - extraxted from dssp files of PDB data base: cullpdb_pc25_res1.6_R0.25_d170113_chains3953.gz
    //Most Frequent ------> Less Frequent

    public final static String[] DsspSTRUCTURESortedByFrequency = {"H_X_S+__","_____-__","H_>_S+__","E____-A_","T3__S+__","E____-AA","_____+__","S___S+__","________","S___S-__","H_<_S+__","E____-a_","_>___-__","E____-aa","E____-_A","E____+A_","H><_S+__","H3>_S+__","_<___-__","T3__S-__","H<>_S+__","H3<_S+__","__>__-__","G>__S+__","G<__S+__","H3<5S+__","E____+AA","H><5S+__","H_<5S+__","G3__S+__","S<__S-__","E____-aA","___<_-__","_<___+__","_>>__-__","E____-_a","T34_S+__","H_<>S+__","E____+_A","H_X>S+__","H>X_S+__","S<__S+__","T3<5S-__","E____+a_","T3<_S+__","E____-__","E____-Aa","T_4_S+__","T<_5S+__","B____-A_","E____+aa","S>__S-__","T<4_S+__","_>___+__","E____+__","T<__S+__","H<<_S+__","__<__-__","E>___-A_","E<___-A_","H3X_S+__","S_>_S-__","__<__+__","S_<_S-__","T<_5_+__","H34_S+__","H<X_S+__","__<_____","H<4_S+__","H_4_S+__","T<<_S+__","H_<_S-__","GX__S+__","G<4_S+__","H_>__+__","H>>_S+__","T3___+__","E___S+__","H_<5S-__","G34_S+__","H><>S+__","__>__+__","G>4_S+__","H_<_____","H3<5S-__","___<_+__","B____-a_","__>_____","E____+_a","T>__S+__","T3>_S+__","T<<5S-__","H_X5S+__","_<______","SX__S-__","E____+aA","T_<5_+__","_X___-__","T_4_S-__","T_45S+__","T>>_S+__","G<>_S+__","E___S-A_","H>4_S+__","S>__S+__","T_<5S+__","E___S-__","T<<5_+__","E>___-_A","S_<_S+__","E<__S-A_","H_><S+__","G>___+__","E>__S-A_","S__<S+__","S>>_S-__","T345S+__","_>>__+__","E____+Aa","T<<5S+__","E___S+A_","S__<S-__","E_>__-A_","T<_5S-__","S<<_S-__","_<>__-__","_><__-__","B>___-A_","G><_S+__","_>______","_X___+__","B____+A_","S_>_S+__","S<>_S-__","_<<__-__","_><__+__","HX4_S+__","T<4_S-__","T>4_S+__","GX___+__","G<_5S-__","E<___+A_","E3__S+__","GX>_S+__","G>>_S+__","E<___-a_","E_____A_","H_>>S+__","E>__S-_A","T_45S-__","E<___-AA","T3_<S+__","G3_5S+__","T3_5S+__","T<_5_-__","HX<_S+__","E>___-a_","E<___-_A","B<___-A_","G<__S-__","B___S-A_","GX4_S+__","HX>_S+__","H<>__+__","E__<_-A_","G>_5S+__","T<___+__","E___S-a_","_<>__+__","SX>_S-__","G3__S-__","E>>__-A_","T3______","G>__S-__","B____+a_","E___S-_A","T<__S-__","T34_S-__","T3<5S+__","E>__S-AA","T3___-__","T__5S+__","S><_S-__","_X>__-__","E>___-_a","E_<__-A_","H_4>S+__","_<<__+__","H_>5S+__","B<__S-A_","T><_S+__","H_45S+__","__><_+__","T_X5S+__","H3>__+__","TX4_S+__","S>>_S+__","B>__S-A_","T<45S-__","_<<_____","G<_5S+__","B_>__-A_","E<___+a_","T3_5S-__","E__<_-AA","GX<_S+__","___<____","T_<5S-__","E______A","E_____AA","__X__-__","_>_<_+__","T<>_S+__","T3<_S-__","G<___+__","HXX_S+__","H>X>S+__","T3X_S+__","H>45S+__","T_<<S+__","G<<_S+__","E_>>_-A_","T_4__+__","E___S-aa","T__5S-__","E___S+_A","E___S+a_","T<4__+__","E>___-AA","H<45S+__","H<<5S+__","H345S+__","B___S-a_","T>X_S+__","G345S+__","S_<<S+__","T>45S+__","T_4_____","_X>__+__","T3<_____","T>__S-__","H<<5S-__","E<__S-_A","G>45S+__","TX>_S+__","B<___+A_","H<>>S+__","__<<_+__","E_>__-_A","E_____a_","G<______","_>_>_-__","G3<_S+__","E>___+A_","E_>__-a_","H<<_S-__","E___S-_a","E>___-__","T><5S-__","B>>__-A_","__X__+__","S><_S+__","T>>__+__","H3>>S+__","T<<_S-__","T<45S+__","T34__+__","SX__S+__","B_____A_","T>___+__","H3<_____","_>>_____","T_>_S+__","B___S+A_","_>_<_-__","_X<__-__","TX<_S+__","_>X__+__","T><5S+__","E>___-aa","E<___-__","H_X<S+__","E___S-AA","E>>__-_A","E<___+AA","H3><S+__","E_>>_-_A","E>>__-a_","__<<_-__","G>_>S+__","E__<S-A_","H<4>S+__","G<45S-__","T345S-__","SX<_S-__","_>X__-__","T3>__+__","E_>>_-AA","B>___-a_","T__<S+__","GX>__+__","E<___+_A","H><_S-__","E_<__+A_","T__5_+__","S<>_S+__","S<<_S+__","T_4>S+__","B_<__-A_","H<<_____","T<4_____","H>4>S+__","E_<__-AA","E_<__-a_","H3>_S-__","T_4<S-__","H_<<S-__","E___S+aa","T_<5_-__","T<<5_-__","G<45S+__","T>_5S+__","G>_<S+__","E>>>_-A_","__><_-__","T<___-__","T<<_____","H<>5S+__","T_4<S+__","E__<_-a_","B__<_-A_","E>___+_A","B<___-a_","TX_5S+__","S_<<S-__","HX45S+__","G<4_S-__","B___S+a_","B<__S-a_","T<>5S+__","G>>__+__","_>>>_-__","TX__S+__","T<X_S+__","G<>__+__","_>><_+__","H>>__+__","E___S+_a","B_____a_","H_>_S-__","H3X>S+__","B_<_S-A_","GX_5S-__","H3>5S+__","E_____aa","GX__S-__","G<<5S-__","E>>>_-_A","T>>_S-__","___>_-__","TX_5_+__","E__<_+A_","T_>5S+__","T>_5S-__","T3_<S-__","T>X5S+__","T__5_-__","T<______","E<__S+A_","E_>__-_a","G><5S+__","H3<_S-__","GX_5S+__","G><__+__","E<__S-AA","E<__S-a_","E>__S+A_","__>>_-__","T<_5____","H34>S+__","E__<_-_A","S<X_S-__","B____-AA","H<<>S+__","E_<__-_A","B_>__-a_","E___S+AA","B345S+a_","_X<__+__","T3<__+__","S_X_S+__","H_<<S+__","E>>__-_a","TX_5S-__","T__<S-__","T_<_S+__","H<X>S+__","G>_<_+__","E<___-Aa","_<X__-__","T3<__-__","HX<5S+__","B<__S+A_","_>_>_+__","TX>__+__","T3_<_+__","G>4__+__","B<___+a_","TXX_S+__","T_X<S+__","T>4_S-__","T34<S+__","S>X_S-__","HX>__+__","E<___+__","T>_<S+__","H3<>S+__","G3_<S+__","E<__S-__","E3__S+A_","H>><S+__","H><5S-__","E>>__-AA","E_____aA","E<___-_a","H_><_+__","E3__S+a_","B3__S+a_","E_<__+a_","H>>_S-__","GX4__+__","E>__S+_A","E__<_-aA","B____+aa","E<___-aA","T3>_S-__","HX<>S+__","HX4>S+__","E3__S-__","T>_5_+__","T<>__+__","G>4>S+__","E_<_S-A_","G>4<S+__","E___S-aA","E__<_+AA","T>4>S+__","S_X_S-__","G<4_____","E__>_-A_","E______a","SXX_S-__","S>_>S-__","G>_5S-__","G3_5S-__","E>>_S-A_","E<<__-A_","B__<S-a_","___>____","TX<5S-__","T_<5____","T<>5_+__","T<4>S+__","T>4__+__","T<<__+__","H3X5S+__","G>_5_+__","E<__S+__","E>_>_-_A","E__<_+a_","__X<_-__","T_><S+__","S>_<S-__","H<X5S+__","H_4<S+__","G>X_S+__","G>___-__","E3__S+aa","B>___+A_","_<X__+__","TX___+__","T<4__-__","H>4<S+__","E<<__+A_","E>___+_a","B_<_S-a_","B__<_+A_","B>>__-a_","_><<_-__","TX>5S+__","G34_S-__","E___S-Aa","E>__S-Aa","E_<_S-a_","B_>_S-A_","__<<____","_XX__-__","T>_>S+__","T><>S+__","GX_>S+__","G><_S-__","E>>>_-AA","E>_>_-A_","B>_>_-A_","B__<_-a_","T>X<S+__","T>><S+__","T34>S+__","H>X5S+__","H_>>_+__","GX<5S+__","E>__S-aA","E__>_-AA","___>_+__","TX<5S+__","T>45S-__","T_4__-__","H>>>S+__","E_>__-AA","H34<S+__","H_>__-__","GX<__+__","E__<_-Aa","B>>_S-A_","B<____A_","B_<__-a_","T_<<S-__","T_4>_+__","H_X_S-__","E_<_S+A_","E>___+a_","BX___-A_","B__<S-A_","B>__S+A_","B3__S-a_","T><<S+__","S>X_S+__","G<X_S+__","GX___-__","EX___-A_","E>___-aA","B_>>_-A_","B_<__+A_","_>>>_+__","_><<_+__","T<X5S+__","T3_>S+__","I__<S+__","I_45S+__","H<>__-__","G><__-__","E___S+aA","B>>>_-A_","B3__S-A_","_XX__+__","TX<5_+__","H_>XS+__","H<>>_+__","GX45S+__","E_____Aa","E<___+Aa","E<___+aA","E_>__+A_","E>_>_+_A","E<<__-a_","B<<__-A_","TX_5_-__","T<<5____","SX>_S+__","SX<_S+__","S_><S-__","I_4<S+__","H_X__-__","E<___+_a","E3__S+_a","BX__S-A_","T>>5S+__","T3<<S+__","S><<S-__","I_4>S+__","G3<5S+__","G3___+__","E>>_S-AA","E__>_+A_","E_>>_+A_","E<<__+a_","B>__S-a_","B_____AA","__>>____","TX>5_+__","T>_5_-__","T><5_+__","T>4<S+__","T3>5S+__","H_<5____","H_45S-__","G><>S+__","E>__S-_a","E>_>_-AA","E_>>_+_A","E_<__+_A","E>_>_+A_","B____+AA","B____-Aa","B3__S+A_","TX4__+__","T3><S+__","T3<<S-__","G>45S-__","G>><_+__","E3__S+_A","B___S-aA","B_>_S-a_","B34_S+a_","__>>_+__","_>><_-__","TX<__+__","S>_>S+__","S>>>S+__","H<4_S-__","GX_5_+__","G>>>S+__","G>4_S-__","G3_<S-__","G345S-__","G>_>_+__","E__<S-AA","E>___+AA","E__>_-_A","E>>__+A_","E>>>_+A_","E__5S+A_","E3__S-a_","T3>>S+__","S__>S-__","GXX_S+__","GX>>S+__","G><<S+__","E>___+__","B345S-a_","TX4_S-__","T_<_S-__","T_<>S+__","T>>>S+__","T_>5S-__","T<>5_-__","T34__-__","T><__+__","SX_>S-__","S>><S-__","I__5S+__","H>X__+__","H><__+__","GX_5_-__","G>>_S-__","G3>_S+__","G34<S+__","E__<S-_A","E>>_S-_A","E>__S-a_","E<<_S-a_","BX__S-a_","B____-aA","B_____aa","B>_>_+A_","TXX__+__","TX<_S-__","T_X>S+__","T>X>S+__","TX45S+__","T__>S+__","T_4XS+__","S_>>S+__","HX>>S+__","H3><_+__","H_<__+__","GX4>S+__","G<>>S+__","E<__S+_A","E_>>_-aA","E__>_+_A","E3__S-A_","E3__S-_a","B____-aa","B__>_+A_","B>4_S-A_","_>_>____","TXX5_+__","T><_S-__","T3_5_+__","T>___-__","I_4XS+__","H>>5S+__","H_4XS+__","H3X<S+__","H3<<S+__","G<<5S+__","G<___-__","G<<_____","E>__S+__","E__>_+AA","E_>>_+AA","E>___-Aa","E__<_-_a","E__5_+A_","E3___+a_","E3___+__","BX___-a_","B__<S-aA","B_<_S+A_","B_<__+a_","_><>_-__","_>X<_-__","TXX5S+__","TX__S-__","T_45_+__","T3<5_-__","T_<<_+__","S<X_S+__","I_<>S+__","H<><S+__","H>45S-__","H<45S-__","H>>__-__","GX>_S-__","GX<>S+__","G<>5S+__","G<4__+__","E___S+Aa","E_<__+AA","E_>__+_a","E_4_S+A_","E3___-__","BX___+A_","B<<_S-a_","B_>__+A_","_X_>_+__","_X_>_-__","TX4>S+__","TX45S-__","T>4>_+__","T345_+__","T34<_+__","T3_<_-__","S>_<S+__","S>>>S-__","I_<5S+__","H>X<S+__","H><<S+__","H_><_-__","GX_>_+__","G>_>S-__","EX___-_A","E_>_S-A_","E>__S-__","E_>>_-Aa","E_<__-Aa","E__<_+_A","E>>__+_A","E_>__+a_","E_<__-_a","B__>_-A_","B_<___A_","B__<_+a_","B>_>_-a_","B3<_S+a_","_<>>_-__","_>X<_+__","TX>>S+__","TX>5S-__","TX>5_-__","T<X5_+__","T>45_+__","T<45_-__","T34<S-__","T34_____","T3<<_-__","I_><S+__","I3_<S+__","H>X_S-__","H34_S-__","H>><_+__","H<<__-__","G>><S+__","G34>S+__","E<__S+AA","E<__S-aA","E_<_S+_A","E_<<_-A_","E__<_+_a","E3__S-_A","BX>_S-A_","BX>__-A_","B<__S-AA","B___S+aa","B<>_S-A_","B<<_S-A_","B_>>_+A_","B>___+a_","B<45S+a_","TX<>S+__","T_X5_+__","T>X__+__","T<X__+__","T>_<S-__","T<>_S-__","T<<>S+__","T_>5_+__","T>>5_-__","T><5_-__","T3_>S-__","T__<____","T>_>_+__","T>_<_+__","T><__-__","S_>>S-__","S>><S+__","S<>>S+__","I__>S+__","HXX>S+__","H_<XS+__","HX>__-__","H<X__+__","H<>5_+__","H_4_S-__","H3X__+__","H345S-__","H3>>_+__","G>4>_+__","G>4<_+__","EX___+A_","E__<S+aa","E_<_S+a_","E_<__-aA","E<___-aa","E__>_-a_","E_<__+_a","BX__S+A_","BX>__-a_","B___S-aa","B><_S-A_","B_<_S+a_","B____+aA","B_45S+a_","_<_>_+__","_<_>_-__","TX>_S-__","T>X5S-__","TX<5_-__","T>_>S-__","T<>>S+__","T>>5S-__","T<>5S-__","T_4<_+__","T3X__+__","T3<>S+__","T34>_+__","T3_<____","T__<_+__","T__<_-__","T_<<_-__","T<<__-__","SX<>S-__","S__>S+__","S_><S+__","S><<S+__","I__5S-__","H_X__+__","H>X__-__","GX_>S-__","GX<5S-__","GX<__-__","G<<_S-__","G>>5S+__","G<4>S+__","G>45_+__","G><<_-__","EX__S-A_","E>__S+AA","E__<S-aA","E_>_S+_A","E_<_S-_A","E>__S+a_","E<__S+a_","E>>__-aA","E>>>_+_A","E<_<_-A_","E__>_-_a","E3___+A_","B>__S+a_","B<__S+a_","B<>_S-a_","B<<_S+a_","B>___-AA","B_____Aa","B><__-A_","B<<__+A_","B<____a_","B<>__-a_","B>_5S-A_","B<_5S-A_","B__5S-a_","B<<5S+a_","B3<_S+A_","B3___+A_","B3<__-A_","B3<5S+a_","_<>>_+__","_X>>_+__","__X<_+__","TX_>S+__","T_X_S+__","T__XS+__","T_<XS+__","T>X_S-__","T><>S-__","T<_>S+__","T>>5_+__","T<45_+__","T>4__-__","T3X>S-__","T3_5_-__","T_<<____","T>>__-__","T>>>_+__","T><>_+__","S_X<S+__","S_<>S-__","S<_>S-__","I_>XS+__","I>><S+__","I><5S+__","I3_5S+__","I34>S+__","I34<S+__","H_XXS+__","H_X>S-__","H_X5S-__","H_X<_-__","H_>>S+A_","H_>>S-__","H_<>S-__","H_4>S-__","H3>>S-__","H_>>_-__","H_<__-__","G>X>S+__","GX4_S-__","G>X__+__","G>_<S-__","G><<S-__","G<_>S+__","G<>_S-__","G<<__+__","EX__S+A_","EX___+_A","E>>__-Aa","E_>__-aA","E_<__-aa","E_>__+_A","E_<<_+A_","E_>>_-_a","E>>>_-a_","E__5S+_A","E_45S+A_","E3__S+AA","E3___-_a","B_X__-A_","B_>_S+A_","B__>S-a_","B__<S+a_","B>>_S-a_","B____+Aa","B__<_-aA","B>___-aA","B_><_+A_","B<>__-A_","B<<___A_","B__<__a_","B<<__-a_","B_4_S+A_","B<4_S+a_","B3<_S-A_","B3_5S+a_","_X>>_-__","_X<>_-__","T>X5_+__","TX4__-__","TX>__-__","T>X__-__","T_>>S+__","T__5____","T>4<S-__","T_45_-__","T3X<S+__","T3><_+__","T3<<_+__","T_<__-__","T>_<_-__","T><<_+__","SXX_S+__","SX>>S-__","S_X<S-__","S><>S-__","S<_>S+__","I_X>S+__","I__XS+__","I_<XS+__","I__<S-__","I_<<S+__","I>_>S+__","I>_<S+__","I<_>S+__","I_>5S+__","I_<5S-__","I>_5S+__","I<_5S+__","I>4XS+__","I<4>S+__","I<45S+__","I3<5S+__","I__>_+__","I__<_-__","HXX__+__","HX><S+__","HX<_S-__","HX<5S-__","HX45S-__","H_>X_+__","H><<S-__","H<>_S-__","H_>5_+__","H>4_S-__","H3X_S-__","H3X>S-__","H3X__-__","H>>>_+__","H><__-__","H<<__+__","GXX__+__","GX>>S-__","GX<_S-__","G>_XS+__","G<X_S-__","GX>__-__","G>_X_+__","G<X__+__","G>4>S-__","G>45_-__","G3>>S+__","G>>__-__","G><>_+__","EX__S-_A","E__<S+AA","E_>_S-AA","E_<_S+AA","E__>S-A_","E_>_S-_A","E>_>S-A_","E__<S-a_","E<__S-_a","E_>__+AA","E>_>_+AA","E__<_+Aa","E__>_-aA","E__<_+aA","E>>>_-aA","E_>__-aa","E>___+aa","E<<__+_A","E<<__-_A","E_>>_-a_","E_<<_-a_","E>_>_-a_","E>>__+a_","E__5_+_A","E__5_-_A","E__5_+__","E_4_S+_A","E_4_S-A_","E3___-AA","E3___+_a","E3_5S+A_","BX_>S-a_","BX___+a_","BX_5S-A_","B___S-AA","B__>S+A_","B_>>S-A_","B_<<S-A_","B>_>S-A_","B>>_S+A_","B>>>S-A_","B__>_+AA","B__<_+AA","B<___-AA","B_<<_+A_","B>_<_-A_","B><__+A_","B_>__+a_","B<>__+a_","B__5S+A_","B<_5S+a_","B__5_+A_","B__5_-A_","B_4_S+a_","B_4__+a_","B3_<S-A_","B3<_S-a_","B3___+a_","B3_5S-a_","B3<5S-a_","__<>_+__","_>>>____","_><>_+__","_<<>_+__","_XX>_+__","__>X_-__","_>_X_-__","TXX5_-__","TXX__-__","TX_>S-__","T_>XS+__","T>_XS+__","T>>XS+__","T><XS+__","T<X_S-__","T<X5S-__","T_X5_-__","T>X5_-__","T<X5_-__","TX___-__","TX<__-__","T>X>_+__","T_>_S-__","T_<>S-__","T>>>S-__","T><<S-__","T<<>S-__","T_>5_-__","T>4XS+__","T>4X_+__","T_4>S-__","T>4>S-__","T<4>S-__","T<4<S+__","T>45_-__","T_4>_-__","T_4<_-__","T_4<____","T>4<_+__","T3X>S+__","T3>XS-__","T3X5S+__","T3X5S-__","T3>>S-__","T3<>S-__","T3>5S-__","T3<5_+__","T34>S-__","T34<_-__","T3_>_+__","T3>>_+__","T>><_+__","T><>_-__","T<_>_+__","SXX>S-__","S_X>S+__","S>X>S-__","S>X<S+__","S>_XS+__","S><>S+__","I>_XS+__","IX<5S+__","I_X5S+__","I>X5S+__","I>_<S-__","I<>5S+__","I_4<S-__","I>4>S+__","I<45S-__","I3_<S-__","I34XS+__","I345S+__","I345S-__","I__<____","I>>>_+__","HXX5S+__","H><XS+__","H<X>S-__","HX>5S+__","H>X5S-__","H_X5_+__","HX4_S-__","HX4<S+__","HX<__+__","H_X>_-__","H>X<_+__","H<X__-__","H_>_S-A_","H>>_S+_a","H_><S-__","H>><S-__","H><>S-__","H<>>S-__","H_>5_-__","H>4XS+__","H<4<S+__","H3>XS+__","H3>_S+A_","H3><S-__","H3>5_+__","H34>S-__","H34<S-__","H3>__-__","H3><_-__","H3<__+__","H_<<____","H<>>_-__","GXX5S+__","G>_XS-__","G>>XS+__","GX>5S+__","G<X5S+__","GX<5_+__","GX45S-__","GX45_-__","GX4__-__","GX4>_+__","GX<>_+__","G>X__-__","G>X>_+__","G>><S-__","G<<>S+__","G><5S-__","G<>5S-__","G>_5_-__","G>>5_+__","G><5_-__","G>4__-__","G3X_S+__","G3_>S+__","G3<_S-__","G34<S-__","G34__+__","G3___-__","G>_>_-__","G>_<_-__","G>>>_+__","G<_>_+__","EX>__-A_","EX___+a_","E_>_S+AA","E_>>S-AA","E_<_S-AA","E<<_S+AA","E__<S-Aa","E_<_S+Aa","E>>_S-Aa","E_<_S+aA","E_<_S-aA","E>>_S-aA","E__<S-aa","E_>_S+A_","E>_>S+A_","E>>_S+A_","E>>>S-A_","E<<_S+A_","E<<_S-_A","E_>_S-a_","E_>_S-_a","E_<_S-_a","E>>_S-a_","E>>>S-a_","E<__S+_a","E_>_S-__","E>>_S+__","E>>>_+AA","E<<__+AA","E<<__-AA","E_<__+Aa","E>>>_-Aa","E<<__-Aa","E__>_+aA","E_<__+aA","E>>>_+aA","E<<__-aA","E__<_-aa","E>_<_+A_","E>_<_-A_","E><__-A_","E>_>_+a_","E>_>_-_a","E>>__+_a","E<>__-a_","E__5S-_A","E__5S+__","E__5_+AA","E__5_-AA","E<_5_-AA","E__5_-A_","E>_5_+A_","E_4_S+a_","E_4_S+_a","E_4_S-a_","E_4__+_a","E_4__-_a","E_45S+_A","E_45S-A_","E3__S+aA","E3__S-aa","E3<_S+a_","E3___-aA","E3___-A_","E3___-a_","E34_S+_a","E34_S-_a","E34__-_a","E345S+AA","E345S+_a","E345S+__","E__>_-__","E__<_-__","E_>__-__","E_>>_-__","E>>__-__","E>>>_+__","BXX_S-a_","BX>>S-A_","BX<_S-A_","B_X_S-A_","BX__S+a_","BX>_S+a_","BX>_S-a_","B_X_S-a_","BX>>_-A_","B>X__-A_","B<X__-A_","BX_5S-a_","BX<5S-a_","BX4_S-A_","B___S+AA","B___S+Aa","B___S-Aa","B>__S-Aa","B>__S-aA","B>_<S-aA","B__<S-aa","B_>>S+A_","B_<<S+A_","B>_>S+A_","B><_S+A_","B_<<S-a_","B>_<S-a_","B><_S+a_","B><_S-a_","B__<_-AA","B_<__-AA","B>_>_-AA","B<____AA","B__<_-Aa","B_>__-Aa","B_<__+Aa","B_<__-Aa","B>_>_-Aa","B<____Aa","B_____aA","B__<_+aA","B__>_-aa","B__<__A_","B>>__+A_","B>>>_+A_","B<>__+A_","B_<<_-a_","B>_>_+a_","B>>__+a_","B><>_-a_","B<<__+a_","B>>5S+aA","B__5S+a_","B_<5S+a_","B_<5S-a_","B>_5S-a_","B<_5S-a_","B>_5_+A_","B>_5_-A_","B<_5_-A_","B__5_+a_","B__5_-a_","B_4_S-A_","B<4_S+A_","B<4_S-A_","B_4_S-a_","B>4_S+a_","B<45S-a_","B3__S+aA","B3>_S-A_","B3<__-AA","B3___-A_","B3____a_","B3_5S+A_","B3<5S-A_","B3_5_-a_","B345_+a_","_<_<_+__"};

        //public final static String[] DsspSTRUCTURESortedByFrequency = {"H__X_S+__", "______-__", "H__>_S+__", "E_____-A_", "T_3__S+__", "E_____-AA", "______+__", "S____S+__", "_________", "S____S-__", "H__<_S+__", "E_____-a_", "__>___-__", "E_____-aa", "E_____-_A", "E_____+A_", "H_><_S+__", "H_3>_S+__", "__<___-__", "T_3__S-__", "H_<>_S+__", "H_3<_S+__", "___>__-__", "G_>__S+__", "G_<__S+__", "H_3<5S+__", "E_____+AA", "H_><5S+__", "H__<5S+__", "G_3__S+__", "S_<__S-__", "E_____-aA", "____<_-__", "__<___+__", "__>>__-__", "E_____-_a", "T_34_S+__", "H__<>S+__", "E_____+_A", "H__X>S+__", "H_>X_S+__", "S_<__S+__", "T_3<5S-__", "E_____+a_", "T_3<_S+__", "E_____-__", "E_____-Aa", "T__4_S+__", "T_<_5S+__", "B_____-A_", "E_____+aa", "S_>__S-__", "T_<4_S+__", "__>___+__", "E_____+__", "T_<__S+__", "H_<<_S+__", "___<__-__", "E_>___-A_", "E_<___-A_", "H_3X_S+__", "S__>_S-__", "___<__+__", "S__<_S-__", "T_<_5_+__", "H_34_S+__", "H_<X_S+__", "___<_____", "H_<4_S+__", "H__4_S+__", "T_<<_S+__", "H__<_S-__", "G_X__S+__", "G_<4_S+__", "H__>__+__", "H_>>_S+__", "T_3___+__", "E____S+__", "H__<5S-__", "G_34_S+__", "H_><>S+__", "___>__+__", "G_>4_S+__", "H__<_____", "H_3<5S-__", "____<_+__", "B_____-a_", "___>_____", "E_____+_a", "T_>__S+__", "T_3>_S+__", "T_<<5S-__", "H__X5S+__", "__<______", "S_X__S-__", "E_____+aA", "T__<5_+__", "__X___-__", "T__4_S-__", "T__45S+__", "T_>>_S+__", "G_<>_S+__", "E____S-A_", "H_>4_S+__", "S_>__S+__", "T__<5S+__", "E____S-__", "T_<<5_+__", "E_>___-_A", "S__<_S+__", "E_<__S-A_", "H__><S+__", "G_>___+__", "E_>__S-A_", "S___<S+__", "S_>>_S-__", "T_345S+__", "__>>__+__", "E_____+Aa", "T_<<5S+__", "E____S+A_", "S___<S-__", "E__>__-A_", "T_<_5S-__", "S_<<_S-__", "__<>__-__", "__><__-__", "B_>___-A_", "G_><_S+__", "__>______", "__X___+__", "B_____+A_", "S__>_S+__", "S_<>_S-__", "__<<__-__", "__><__+__", "H_X4_S+__", "T_<4_S-__", "T_>4_S+__", "G_X___+__", "G_<_5S-__", "E_<___+A_", "E_3__S+__", "G_X>_S+__", "G_>>_S+__", "E_<___-a_", "E______A_", "H__>>S+__", "E_>__S-_A", "T__45S-__", "E_<___-AA", "T_3_<S+__", "G_3_5S+__", "T_3_5S+__", "T_<_5_-__", "H_X<_S+__", "E_>___-a_", "E_<___-_A", "B_<___-A_", "G_<__S-__", "B____S-A_", "G_X4_S+__", "H_X>_S+__", "H_<>__+__", "E___<_-A_", "G_>_5S+__", "T_<___+__", "E____S-a_", "__<>__+__", "S_X>_S-__", "G_3__S-__", "E_>>__-A_", "T_3______", "G_>__S-__", "B_____+a_", "E____S-_A", "T_<__S-__", "T_34_S-__", "T_3<5S+__", "E_>__S-AA", "T_3___-__", "T___5S+__", "S_><_S-__", "__X>__-__", "E_>___-_a", "E__<__-A_", "H__4>S+__", "__<<__+__", "H__>5S+__", "B_<__S-A_", "T_><_S+__", "H__45S+__", "___><_+__", "T__X5S+__", "H_3>__+__", "T_X4_S+__", "S_>>_S+__", "B_>__S-A_", "T_<45S-__", "__<<_____", "G_<_5S+__", "B__>__-A_", "E_<___+a_", "T_3_5S-__", "E___<_-AA", "G_X<_S+__", "____<____", "T__<5S-__", "E_______A", "E______AA", "___X__-__", "__>_<_+__", "T_<>_S+__", "T_3<_S-__", "G_<___+__", "H_XX_S+__", "H_>X>S+__", "T_3X_S+__", "H_>45S+__", "T__<<S+__", "G_<<_S+__", "E__>>_-A_", "T__4__+__", "E____S-aa", "T___5S-__", "E____S+_A", "E____S+a_", "T_<4__+__", "E_>___-AA", "H_<45S+__", "H_<<5S+__", "H_345S+__", "B____S-a_", "T_>X_S+__", "G_345S+__", "S__<<S+__", "T_>45S+__", "T__4_____", "__X>__+__", "T_3<_____", "T_>__S-__", "H_<<5S-__", "E_<__S-_A", "G_>45S+__", "T_X>_S+__", "B_<___+A_", "H_<>>S+__", "___<<_+__", "E__>__-_A", "E______a_", "G_<______", "__>_>_-__", "G_3<_S+__", "E_>___+A_", "E__>__-a_", "H_<<_S-__", "E____S-_a", "E_>___-__", "T_><5S-__", "B_>>__-A_", "___X__+__", "S_><_S+__", "T_>>__+__", "H_3>>S+__", "T_<<_S-__", "T_<45S+__", "T_34__+__", "S_X__S+__", "B______A_", "T_>___+__", "H_3<_____", "__>>_____", "T__>_S+__", "B____S+A_", "__>_<_-__", "__X<__-__", "T_X<_S+__", "__>X__+__", "T_><5S+__", "E_>___-aa", "E_<___-__", "H__X<S+__", "E____S-AA", "E_>>__-_A", "E_<___+AA", "H_3><S+__", "E__>>_-_A", "E_>>__-a_", "___<<_-__", "G_>_>S+__", "E___<S-A_", "H_<4>S+__", "G_<45S-__", "T_345S-__", "S_X<_S-__", "__>X__-__", "T_3>__+__", "E__>>_-AA", "B_>___-a_", "T___<S+__", "G_X>__+__", "E_<___+_A", "H_><_S-__", "E__<__+A_", "T___5_+__", "S_<>_S+__", "S_<<_S+__", "T__4>S+__", "B__<__-A_", "H_<<_____", "T_<4_____", "H_>4>S+__", "E__<__-AA", "E__<__-a_", "H_3>_S-__", "T__4<S-__", "H__<<S-__", "E____S+aa", "T__<5_-__", "T_<<5_-__", "G_<45S+__", "T_>_5S+__", "G_>_<S+__", "E_>>>_-A_", "___><_-__", "T_<___-__", "T_<<_____", "H_<>5S+__", "T__4<S+__", "E___<_-a_", "B___<_-A_", "E_>___+_A", "B_<___-a_", "T_X_5S+__", "S__<<S-__", "H_X45S+__", "G_<4_S-__", "B____S+a_", "B_<__S-a_", "T_<>5S+__", "G_>>__+__", "__>>>_-__", "T_X__S+__", "T_<X_S+__", "G_<>__+__", "__>><_+__", "H_>>__+__", "E____S+_a", "B______a_", "H__>_S-__", "H_3X>S+__", "B__<_S-A_", "G_X_5S-__", "H_3>5S+__", "E______aa", "G_X__S-__", "G_<<5S-__", "E_>>>_-_A", "T_>>_S-__", "____>_-__", "T_X_5_+__", "E___<_+A_", "T__>5S+__", "T_>_5S-__", "T_3_<S-__", "T_>X5S+__", "T___5_-__", "T_<______", "E_<__S+A_", "E__>__-_a", "G_><5S+__", "H_3<_S-__", "G_X_5S+__", "G_><__+__", "E_<__S-AA", "E_<__S-a_", "E_>__S+A_", "___>>_-__", "T_<_5____", "H_34>S+__", "E___<_-_A", "S_<X_S-__", "B_____-AA", "H_<<>S+__", "E__<__-_A", "B__>__-a_", "E____S+AA", "B_345S+a_", "__X<__+__", "T_3<__+__", "S__X_S+__", "H__<<S+__", "E_>>__-_a", "T_X_5S-__", "T___<S-__", "T__<_S+__", "H_<X>S+__", "G_>_<_+__", "E_<___-Aa", "__<X__-__", "T_3<__-__", "H_X<5S+__", "B_<__S+A_", "__>_>_+__", "T_X>__+__", "T_3_<_+__", "G_>4__+__", "B_<___+a_", "T_XX_S+__", "T__X<S+__", "T_>4_S-__", "T_34<S+__", "S_>X_S-__", "H_X>__+__", "E_<___+__", "T_>_<S+__", "H_3<>S+__", "G_3_<S+__", "E_<__S-__", "E_3__S+A_", "H_>><S+__", "H_><5S-__", "E_>>__-AA", "E______aA", "E_<___-_a", "H__><_+__", "E_3__S+a_", "B_3__S+a_", "E__<__+a_", "H_>>_S-__", "G_X4__+__", "E_>__S+_A", "E___<_-aA", "B_____+aa", "E_<___-aA", "T_3>_S-__", "H_X<>S+__", "H_X4>S+__", "E_3__S-__", "T_>_5_+__", "T_<>__+__", "G_>4>S+__", "E__<_S-A_", "G_>4<S+__", "E____S-aA", "E___<_+AA", "T_>4>S+__", "S__X_S-__", "G_<4_____", "E___>_-A_", "E_______a", "S_XX_S-__", "S_>_>S-__", "G_>_5S-__", "G_3_5S-__", "E_>>_S-A_", "E_<<__-A_", "B___<S-a_", "____>____", "T_X<5S-__", "T__<5____", "T_<>5_+__", "T_<4>S+__", "T_>4__+__", "T_<<__+__", "H_3X5S+__", "G_>_5_+__", "E_<__S+__", "E_>_>_-_A", "E___<_+a_", "___X<_-__", "T__><S+__", "S_>_<S-__", "H_<X5S+__", "H__4<S+__", "G_>X_S+__", "G_>___-__", "E_3__S+aa", "B_>___+A_", "__<X__+__", "T_X___+__", "T_<4__-__", "H_>4<S+__", "E_<<__+A_", "E_>___+_a", "B__<_S-a_", "B___<_+A_", "B_>>__-a_", "__><<_-__", "T_X>5S+__", "G_34_S-__", "E____S-Aa", "E_>__S-Aa", "E__<_S-a_", "B__>_S-A_", "___<<____", "__XX__-__", "T_>_>S+__", "T_><>S+__", "G_X_>S+__", "G_><_S-__", "E_>>>_-AA", "E_>_>_-A_", "B_>_>_-A_", "B___<_-a_", "T_>X<S+__", "T_>><S+__", "T_34>S+__", "H_>X5S+__", "H__>>_+__", "G_X<5S+__", "E_>__S-aA", "E___>_-AA", "____>_+__", "T_X<5S+__", "T_>45S-__", "T__4__-__", "H_>>>S+__", "E__>__-AA", "H_34<S+__", "H__>__-__", "G_X<__+__", "E___<_-Aa", "B_>>_S-A_", "B_<____A_", "B__<__-a_", "T__<<S-__", "T__4>_+__", "H__X_S-__", "E__<_S+A_", "E_>___+a_", "B_X___-A_", "B___<S-A_", "B_>__S+A_", "B_3__S-a_", "T_><<S+__", "S_>X_S+__", "G_<X_S+__", "G_X___-__", "E_X___-A_", "E_>___-aA", "B__>>_-A_", "B__<__+A_", "__>>>_+__", "__><<_+__", "T_<X5S+__", "T_3_>S+__", "I___<S+__", "I__45S+__", "H_<>__-__", "G_><__-__", "E____S+aA", "B_>>>_-A_", "B_3__S-A_", "__XX__+__", "T_X<5_+__", "H__>XS+__", "H_<>>_+__", "G_X45S+__", "E______Aa", "E_<___+Aa", "E_<___+aA", "E__>__+A_", "E_>_>_+_A", "E_<<__-a_", "B_<<__-A_", "T_X_5_-__", "T_<<5____", "S_X>_S+__", "S_X<_S+__", "S__><S-__", "I__4<S+__", "H__X__-__", "E_<___+_a", "E_3__S+_a", "B_X__S-A_", "T_>>5S+__", "T_3<<S+__", "S_><<S-__", "I__4>S+__", "G_3<5S+__", "G_3___+__", "E_>>_S-AA", "E___>_+A_", "E__>>_+A_", "E_<<__+a_", "B_>__S-a_", "B______AA", "___>>____", "T_X>5_+__", "T_>_5_-__", "T_><5_+__", "T_>4<S+__", "T_3>5S+__", "H__<5____", "H__45S-__", "G_><>S+__", "E_>__S-_a", "E_>_>_-AA", "E__>>_+_A", "E__<__+_A", "E_>_>_+A_", "B_____+AA", "B_____-Aa", "B_3__S+A_", "T_X4__+__", "T_3><S+__", "T_3<<S-__", "G_>45S-__", "G_>><_+__", "E_3__S+_A", "B____S-aA", "B__>_S-a_", "B_34_S+a_", "___>>_+__", "__>><_-__", "T_X<__+__", "S_>_>S+__", "S_>>>S+__", "H_<4_S-__", "G_X_5_+__", "G_>>>S+__", "G_>4_S-__", "G_3_<S-__", "G_345S-__", "G_>_>_+__", "E___<S-AA", "E_>___+AA", "E___>_-_A", "E_>>__+A_", "E_>>>_+A_", "E___5S+A_", "E_3__S-a_", "T_3>>S+__", "S___>S-__", "G_XX_S+__", "G_X>>S+__", "G_><<S+__", "E_>___+__", "B_345S-a_", "T_X4_S-__", "T__<_S-__", "T__<>S+__", "T_>>>S+__", "T__>5S-__", "T_<>5_-__", "T_34__-__", "T_><__+__", "S_X_>S-__", "S_>><S-__", "I___5S+__", "H_>X__+__", "H_><__+__", "G_X_5_-__", "G_>>_S-__", "G_3>_S+__", "G_34<S+__", "E___<S-_A", "E_>>_S-_A", "E_>__S-a_", "E_<<_S-a_", "B_X__S-a_", "B_____-aA", "B______aa", "B_>_>_+A_", "T_XX__+__", "T_X<_S-__", "T__X>S+__", "T_>X>S+__", "T_X45S+__", "T___>S+__", "T__4XS+__", "S__>>S+__", "H_X>>S+__", "H_3><_+__", "H__<__+__", "G_X4>S+__", "G_<>>S+__", "E_<__S+_A", "E__>>_-aA", "E___>_+_A", "E_3__S-A_", "E_3__S-_a", "B_____-aa", "B___>_+A_", "B_>4_S-A_", "__>_>____", "T_XX5_+__", "T_><_S-__", "T_3_5_+__", "T_>___-__", "I__4XS+__", "H_>>5S+__", "H__4XS+__", "H_3X<S+__", "H_3<<S+__", "G_<<5S+__", "G_<___-__", "G_<<_____", "E_>__S+__", "E___>_+AA", "E__>>_+AA", "E_>___-Aa", "E___<_-_a", "E___5_+A_", "E_3___+a_", "E_3___+__", "B_X___-a_", "B___<S-aA", "B__<_S+A_", "B__<__+a_", "__><>_-__", "__>X<_-__", "T_XX5S+__", "T_X__S-__", "T__45_+__", "T_3<5_-__", "T__<<_+__", "S_<X_S+__", "I__<>S+__", "H_<><S+__", "H_>45S-__", "H_<45S-__", "H_>>__-__", "G_X>_S-__", "G_X<>S+__", "G_<>5S+__", "G_<4__+__", "E____S+Aa", "E__<__+AA", "E__>__+_a", "E__4_S+A_", "E_3___-__", "B_X___+A_", "B_<<_S-a_", "B__>__+A_", "__X_>_+__", "__X_>_-__", "T_X4>S+__", "T_X45S-__", "T_>4>_+__", "T_345_+__", "T_34<_+__", "T_3_<_-__", "S_>_<S+__", "S_>>>S-__", "I__<5S+__", "H_>X<S+__", "H_><<S+__", "H__><_-__", "G_X_>_+__", "G_>_>S-__", "E_X___-_A", "E__>_S-A_", "E_>__S-__", "E__>>_-Aa", "E__<__-Aa", "E___<_+_A", "E_>>__+_A", "E__>__+a_", "E__<__-_a", "B___>_-A_", "B__<___A_", "B___<_+a_", "B_>_>_-a_", "B_3<_S+a_", "__<>>_-__", "__>X<_+__", "T_X>>S+__", "T_X>5S-__", "T_X>5_-__", "T_<X5_+__", "T_>45_+__", "T_<45_-__", "T_34<S-__", "T_34_____", "T_3<<_-__", "I__><S+__", "I_3_<S+__", "H_>X_S-__", "H_34_S-__", "H_>><_+__", "H_<<__-__", "G_>><S+__", "G_34>S+__", "E_<__S+AA", "E_<__S-aA", "E__<_S+_A", "E__<<_-A_", "E___<_+_a", "E_3__S-_A", "B_X>_S-A_", "B_X>__-A_", "B_<__S-AA", "B____S+aa", "B_<>_S-A_", "B_<<_S-A_", "B__>>_+A_", "B_>___+a_", "B_<45S+a_", "T_X<>S+__", "T__X5_+__", "T_>X__+__", "T_<X__+__", "T_>_<S-__", "T_<>_S-__", "T_<<>S+__", "T__>5_+__", "T_>>5_-__", "T_><5_-__", "T_3_>S-__", "T___<____", "T_>_>_+__", "T_>_<_+__", "T_><__-__", "S__>>S-__", "S_>><S+__", "S_<>>S+__", "I___>S+__", "H_XX>S+__", "H__<XS+__", "H_X>__-__", "H_<X__+__", "H_<>5_+__", "H__4_S-__", "H_3X__+__", "H_345S-__", "H_3>>_+__", "G_>4>_+__", "G_>4<_+__", "E_X___+A_", "E___<S+aa", "E__<_S+a_", "E__<__-aA", "E_<___-aa", "E___>_-a_", "E__<__+_a", "B_X__S+A_", "B_X>__-a_", "B____S-aa", "B_><_S-A_", "B__<_S+a_", "B_____+aA", "B__45S+a_", "__<_>_+__", "__<_>_-__", "T_X>_S-__", "T_>X5S-__", "T_X<5_-__", "T_>_>S-__", "T_<>>S+__", "T_>>5S-__", "T_<>5S-__", "T__4<_+__", "T_3X__+__", "T_3<>S+__", "T_34>_+__", "T_3_<____", "T___<_+__", "T___<_-__", "T__<<_-__", "T_<<__-__", "S_X<>S-__", "S___>S+__", "S__><S+__", "S_><<S+__", "I___5S-__", "H__X__+__", "H_>X__-__", "G_X_>S-__", "G_X<5S-__", "G_X<__-__", "G_<<_S-__", "G_>>5S+__", "G_<4>S+__", "G_>45_+__", "G_><<_-__", "E_X__S-A_", "E_>__S+AA", "E___<S-aA", "E__>_S+_A", "E__<_S-_A", "E_>__S+a_", "E_<__S+a_", "E_>>__-aA", "E_>>>_+_A", "E_<_<_-A_", "E___>_-_a", "E_3___+A_", "B_>__S+a_", "B_<__S+a_", "B_<>_S-a_", "B_<<_S+a_", "B_>___-AA", "B______Aa", "B_><__-A_", "B_<<__+A_", "B_<____a_", "B_<>__-a_", "B_>_5S-A_", "B_<_5S-A_", "B___5S-a_", "B_<<5S+a_", "B_3<_S+A_", "B_3___+A_", "B_3<__-A_", "B_3<5S+a_", "__<>>_+__", "__X>>_+__", "___X<_+__", "T_X_>S+__", "T__X_S+__", "T___XS+__", "T__<XS+__", "T_>X_S-__", "T_><>S-__", "T_<_>S+__", "T_>>5_+__", "T_<45_+__", "T_>4__-__", "T_3X>S-__", "T_3_5_-__", "T__<<____", "T_>>__-__", "T_>>>_+__", "T_><>_+__", "S__X<S+__", "S__<>S-__", "S_<_>S-__", "I__>XS+__", "I_>><S+__", "I_><5S+__", "I_3_5S+__", "I_34>S+__", "I_34<S+__", "H__XXS+__", "H__X>S-__", "H__X5S-__", "H__X<_-__", "H__>>S+A_", "H__>>S-__", "H__<>S-__", "H__4>S-__", "H_3>>S-__", "H__>>_-__", "H__<__-__", "G_>X>S+__", "G_X4_S-__", "G_>X__+__", "G_>_<S-__", "G_><<S-__", "G_<_>S+__", "G_<>_S-__", "G_<<__+__", "E_X__S+A_", "E_X___+_A", "E_>>__-Aa", "E__>__-aA", "E__<__-aa", "E__>__+_A", "E__<<_+A_", "E__>>_-_a", "E_>>>_-a_", "E___5S+_A", "E__45S+A_", "E_3__S+AA", "E_3___-_a", "B__X__-A_", "B__>_S+A_", "B___>S-a_", "B___<S+a_", "B_>>_S-a_", "B_____+Aa", "B___<_-aA", "B_>___-aA", "B__><_+A_", "B_<>__-A_", "B_<<___A_", "B___<__a_", "B_<<__-a_", "B__4_S+A_", "B_<4_S+a_", "B_3<_S-A_", "B_3_5S+a_", "__X>>_-__", "__X<>_-__", "T_>X5_+__", "T_X4__-__", "T_X>__-__", "T_>X__-__", "T__>>S+__", "T___5____", "T_>4<S-__", "T__45_-__", "T_3X<S+__", "T_3><_+__", "T_3<<_+__", "T__<__-__", "T_>_<_-__", "T_><<_+__", "S_XX_S+__", "S_X>>S-__", "S__X<S-__", "S_><>S-__", "S_<_>S+__", "I__X>S+__", "I___XS+__", "I__<XS+__", "I___<S-__", "I__<<S+__", "I_>_>S+__", "I_>_<S+__", "I_<_>S+__", "I__>5S+__", "I__<5S-__", "I_>_5S+__", "I_<_5S+__", "I_>4XS+__", "I_<4>S+__", "I_<45S+__", "I_3<5S+__", "I___>_+__", "I___<_-__", "H_XX__+__", "H_X><S+__", "H_X<_S-__", "H_X<5S-__", "H_X45S-__", "H__>X_+__", "H_><<S-__", "H_<>_S-__", "H__>5_+__", "H_>4_S-__", "H_3X_S-__", "H_3X>S-__", "H_3X__-__", "H_>>>_+__", "H_><__-__", "H_<<__+__", "G_XX__+__", "G_X>>S-__", "G_X<_S-__", "G_>_XS+__", "G_<X_S-__", "G_X>__-__", "G_>_X_+__", "G_<X__+__", "G_>4>S-__", "G_>45_-__", "G_3>>S+__", "G_>>__-__", "G_><>_+__", "E_X__S-_A", "E___<S+AA", "E__>_S-AA", "E__<_S+AA", "E___>S-A_", "E__>_S-_A", "E_>_>S-A_", "E___<S-a_", "E_<__S-_a", "E__>__+AA", "E_>_>_+AA", "E___<_+Aa", "E___>_-aA", "E___<_+aA", "E_>>>_-aA", "E__>__-aa", "E_>___+aa", "E_<<__+_A", "E_<<__-_A", "E__>>_-a_", "E__<<_-a_", "E_>_>_-a_", "E_>>__+a_", "E___5_+_A", "E___5_-_A", "E___5_+__", "E__4_S+_A", "E__4_S-A_", "E_3___-AA", "E_3___+_a", "E_3_5S+A_", "B_X_>S-a_", "B_X___+a_", "B_X_5S-A_", "B____S-AA", "B___>S+A_", "B__>>S-A_", "B__<<S-A_", "B_>_>S-A_", "B_>>_S+A_", "B_>>>S-A_", "B___>_+AA", "B___<_+AA", "B_<___-AA", "B__<<_+A_", "B_>_<_-A_", "B_><__+A_", "B__>__+a_", "B_<>__+a_", "B___5S+A_", "B_<_5S+a_", "B___5_+A_", "B___5_-A_", "B__4_S+a_", "B__4__+a_", "B_3_<S-A_", "B_3<_S-a_", "B_3___+a_", "B_3_5S-a_", "B_3<5S-a_", "___<>_+__", "__>>>____", "__><>_+__", "__<<>_+__", "__XX>_+__", "___>X_-__", "__>_X_-__", "T_XX5_-__", "T_XX__-__", "T_X_>S-__", "T__>XS+__", "T_>_XS+__", "T_>>XS+__", "T_><XS+__", "T_<X_S-__", "T_<X5S-__", "T__X5_-__", "T_>X5_-__", "T_<X5_-__", "T_X___-__", "T_X<__-__", "T_>X>_+__", "T__>_S-__", "T__<>S-__", "T_>>>S-__", "T_><<S-__", "T_<<>S-__", "T__>5_-__", "T_>4XS+__", "T_>4X_+__", "T__4>S-__", "T_>4>S-__", "T_<4>S-__", "T_<4<S+__", "T_>45_-__", "T__4>_-__", "T__4<_-__", "T__4<____", "T_>4<_+__", "T_3X>S+__", "T_3>XS-__", "T_3X5S+__", "T_3X5S-__", "T_3>>S-__", "T_3<>S-__", "T_3>5S-__", "T_3<5_+__", "T_34>S-__", "T_34<_-__", "T_3_>_+__", "T_3>>_+__", "T_>><_+__", "T_><>_-__", "T_<_>_+__", "S_XX>S-__", "S__X>S+__", "S_>X>S-__", "S_>X<S+__", "S_>_XS+__", "S_><>S+__", "I_>_XS+__", "I_X<5S+__", "I__X5S+__", "I_>X5S+__", "I_>_<S-__", "I_<>5S+__", "I__4<S-__", "I_>4>S+__", "I_<45S-__", "I_3_<S-__", "I_34XS+__", "I_345S+__", "I_345S-__", "I___<____", "I_>>>_+__", "H_XX5S+__", "H_><XS+__", "H_<X>S-__", "H_X>5S+__", "H_>X5S-__", "H__X5_+__", "H_X4_S-__", "H_X4<S+__", "H_X<__+__", "H__X>_-__", "H_>X<_+__", "H_<X__-__", "H__>_S-A_", "H_>>_S+_a", "H__><S-__", "H_>><S-__", "H_><>S-__", "H_<>>S-__", "H__>5_-__", "H_>4XS+__", "H_<4<S+__", "H_3>XS+__", "H_3>_S+A_", "H_3><S-__", "H_3>5_+__", "H_34>S-__", "H_34<S-__", "H_3>__-__", "H_3><_-__", "H_3<__+__", "H__<<____", "H_<>>_-__", "G_XX5S+__", "G_>_XS-__", "G_>>XS+__", "G_X>5S+__", "G_<X5S+__", "G_X<5_+__", "G_X45S-__", "G_X45_-__", "G_X4__-__", "G_X4>_+__", "G_X<>_+__", "G_>X__-__", "G_>X>_+__", "G_>><S-__", "G_<<>S+__", "G_><5S-__", "G_<>5S-__", "G_>_5_-__", "G_>>5_+__", "G_><5_-__", "G_>4__-__", "G_3X_S+__", "G_3_>S+__", "G_3<_S-__", "G_34<S-__", "G_34__+__", "G_3___-__", "G_>_>_-__", "G_>_<_-__", "G_>>>_+__", "G_<_>_+__", "E_X>__-A_", "E_X___+a_", "E__>_S+AA", "E__>>S-AA", "E__<_S-AA", "E_<<_S+AA", "E___<S-Aa", "E__<_S+Aa", "E_>>_S-Aa", "E__<_S+aA", "E__<_S-aA", "E_>>_S-aA", "E___<S-aa", "E__>_S+A_", "E_>_>S+A_", "E_>>_S+A_", "E_>>>S-A_", "E_<<_S+A_", "E_<<_S-_A", "E__>_S-a_", "E__>_S-_a", "E__<_S-_a", "E_>>_S-a_", "E_>>>S-a_", "E_<__S+_a", "E__>_S-__", "E_>>_S+__", "E_>>>_+AA", "E_<<__+AA", "E_<<__-AA", "E__<__+Aa", "E_>>>_-Aa", "E_<<__-Aa", "E___>_+aA", "E__<__+aA", "E_>>>_+aA", "E_<<__-aA", "E___<_-aa", "E_>_<_+A_", "E_>_<_-A_", "E_><__-A_", "E_>_>_+a_", "E_>_>_-_a", "E_>>__+_a", "E_<>__-a_", "E___5S-_A", "E___5S+__", "E___5_+AA", "E___5_-AA", "E_<_5_-AA", "E___5_-A_", "E_>_5_+A_", "E__4_S+a_", "E__4_S+_a", "E__4_S-a_", "E__4__+_a", "E__4__-_a", "E__45S+_A", "E__45S-A_", "E_3__S+aA", "E_3__S-aa", "E_3<_S+a_", "E_3___-aA", "E_3___-A_", "E_3___-a_", "E_34_S+_a", "E_34_S-_a", "E_34__-_a", "E_345S+AA", "E_345S+_a", "E_345S+__", "E___>_-__", "E___<_-__", "E__>__-__", "E__>>_-__", "E_>>__-__", "E_>>>_+__", "B_XX_S-a_", "B_X>>S-A_", "B_X<_S-A_", "B__X_S-A_", "B_X__S+a_", "B_X>_S+a_", "B_X>_S-a_", "B__X_S-a_", "B_X>>_-A_", "B_>X__-A_", "B_<X__-A_", "B_X_5S-a_", "B_X<5S-a_", "B_X4_S-A_", "B____S+AA", "B____S+Aa", "B____S-Aa", "B_>__S-Aa", "B_>__S-aA", "B_>_<S-aA", "B___<S-aa", "B__>>S+A_", "B__<<S+A_", "B_>_>S+A_", "B_><_S+A_", "B__<<S-a_", "B_>_<S-a_", "B_><_S+a_", "B_><_S-a_", "B___<_-AA", "B__<__-AA", "B_>_>_-AA", "B_<____AA", "B___<_-Aa", "B__>__-Aa", "B__<__+Aa", "B__<__-Aa", "B_>_>_-Aa", "B_<____Aa", "B______aA", "B___<_+aA", "B___>_-aa", "B___<__A_", "B_>>__+A_", "B_>>>_+A_", "B_<>__+A_", "B__<<_-a_", "B_>_>_+a_", "B_>>__+a_", "B_><>_-a_", "B_<<__+a_", "B_>>5S+aA", "B___5S+a_", "B__<5S+a_", "B__<5S-a_", "B_>_5S-a_", "B_<_5S-a_", "B_>_5_+A_", "B_>_5_-A_", "B_<_5_-A_", "B___5_+a_", "B___5_-a_", "B__4_S-A_", "B_<4_S+A_", "B_<4_S-A_", "B__4_S-a_", "B_>4_S+a_", "B_<45S-a_", "B_3__S+aA", "B_3>_S-A_", "B_3<__-AA", "B_3___-A_", "B_3____a_", "B_3_5S+A_", "B_3<5S-A_", "B_3_5_-a_", "B_345_+a_", "__<_<_+__"};

    /*private final static LocalStructure[] Dssp30MostFrequentStructures = { new LocalStructure("G3__S+__"),
            new LocalStructure("H_<5S+__"),
            new LocalStructure("H><5S+__"),
            new LocalStructure("E____+AA"),
            new LocalStructure("H3<5S+__"),
            new LocalStructure("G<__S+__"),
            new LocalStructure("G>__S+__"),
            new LocalStructure("__>__-__"),
            new LocalStructure("H3<_S+__"),
            new LocalStructure("H<>_S+__"),
            new LocalStructure("T3__S-__"),
            new LocalStructure("_<___-__"),
            new LocalStructure("H3>_S+__"),
            new LocalStructure("H><_S+__"),
            new LocalStructure("E____+A_"),
            new LocalStructure("E____-_A"),
            new LocalStructure("E____-aa"),
            new LocalStructure("_>___-__"),
            new LocalStructure("E____-a_"),
            new LocalStructure("H_<_S+__"),
            new LocalStructure("S___S-__"),
            new LocalStructure("________"),
            new LocalStructure("S___S+__"),
            new LocalStructure("_____+__"),
            new LocalStructure("E____-AA"),
            new LocalStructure("T3__S+__"),
            new LocalStructure("E____-A_"),
            new LocalStructure("H_>_S+__"),
            new LocalStructure("_____-__"),
            new LocalStructure("H_X_S+__")};

    //private final static LocalStructure[] Dssp100MostFrequentStructures = {};
    private final static LocalStructure[] Dssp100MostFrequentStructures = {new LocalStructure("H><_S+__"),
            new LocalStructure("T_<5_+__"),
            new LocalStructure("E____+aA"),
            new LocalStructure("SX__S-__"),
            new LocalStructure("_<______"),
            new LocalStructure("H_X_S+__"),
            new LocalStructure("H_X5S+__"),
            new LocalStructure("T<<5S-__"),
            new LocalStructure("T3>_S+__"),
            new LocalStructure("T>__S+__"),
            new LocalStructure("E____+_a"),
            new LocalStructure("__>_____"),
            new LocalStructure("E____+A_"),
            new LocalStructure("B____-a_"),
            new LocalStructure("___<_+__"),
            new LocalStructure("H3<5S-__"),
            new LocalStructure("H_<_____"),
            new LocalStructure("G>4_S+__"),
            new LocalStructure("E____-_A"),
            new LocalStructure("E____-aa"),
            new LocalStructure("__>__+__"),
            new LocalStructure("H><>S+__"),
            new LocalStructure("G34_S+__"),
            new LocalStructure("_>___-__"),
            new LocalStructure("H_<5S-__"),
            new LocalStructure("E___S+__"),
            new LocalStructure("T3___+__"),
            new LocalStructure("H>>_S+__"),
            new LocalStructure("H_>__+__"),
            new LocalStructure("G<4_S+__"),
            new LocalStructure("GX__S+__"),
            new LocalStructure("H_<_S-__"),
            new LocalStructure("E____-a_"),
            new LocalStructure("T<<_S+__"),
            new LocalStructure("H_4_S+__"),
            new LocalStructure("H<4_S+__"),
            new LocalStructure("__<_____"),
            new LocalStructure("H<X_S+__"),
            new LocalStructure("H34_S+__"),
            new LocalStructure("H_<_S+__"),
            new LocalStructure("T<_5_+__"),
            new LocalStructure("S_<_S-__"),
            new LocalStructure("__<__+__"),
            new LocalStructure("S_>_S-__"),
            new LocalStructure("H3X_S+__"),
            new LocalStructure("E<___-A_"),
            new LocalStructure("S___S-__"),
            new LocalStructure("________"),
            new LocalStructure("E>___-A_"),
            new LocalStructure("__<__-__"),
            new LocalStructure("H<<_S+__"),
            new LocalStructure("T<__S+__"),
            new LocalStructure("E____+__"),
            new LocalStructure("_>___+__"),
            new LocalStructure("T<4_S+__"),
            new LocalStructure("S___S+__"),
            new LocalStructure("S>__S-__"),
            new LocalStructure("_____+__"),
            new LocalStructure("E____+aa"),
            new LocalStructure("B____-A_"),
            new LocalStructure("T<_5S+__"),
            new LocalStructure("T_4_S+__"),
            new LocalStructure("E____-AA"),
            new LocalStructure("E____-Aa"),
            new LocalStructure("E____-__"),
            new LocalStructure("T3<_S+__"),
            new LocalStructure("E____+a_"),
            new LocalStructure("T3<5S-__"),
            new LocalStructure("S<__S+__"),
            new LocalStructure("H>X_S+__"),
            new LocalStructure("H_X>S+__"),
            new LocalStructure("E____+_A"),
            new LocalStructure("H_<>S+__"),
            new LocalStructure("T3__S+__"),
            new LocalStructure("T34_S+__"),
            new LocalStructure("E____-_a"),
            new LocalStructure("_>>__-__"),
            new LocalStructure("E____-A_"),
            new LocalStructure("_<___+__"),
            new LocalStructure("H_>_S+__"),
            new LocalStructure("___<_-__"),
            new LocalStructure("E____-aA"),
            new LocalStructure("S<__S-__"),
            new LocalStructure("G3__S+__"),
            new LocalStructure("H_<5S+__"),
            new LocalStructure("H><5S+__"),
            new LocalStructure("E____+AA"),
            new LocalStructure("H3<5S+__"),
            new LocalStructure("G<__S+__"),
            new LocalStructure("G>__S+__"),
            new LocalStructure("_____-__"),
            new LocalStructure("__>__-__"),
            new LocalStructure("H3<_S+__"),
            new LocalStructure("H<>_S+__"),
            new LocalStructure("T3__S-__"),
            new LocalStructure("_<___-__"),
            new LocalStructure("H3>_S+__"),
            new LocalStructure("T_45S+__"),
            new LocalStructure("T_4_S-__"),
            new LocalStructure("_X___-__")};
*/

}
