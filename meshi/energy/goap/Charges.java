package meshi.energy.goap;

/**
 * Created by chen on 05/07/2015.
 */
public class Charges {
        public String[][] cind_atomsByOrder;
        float[][] chg_chargesByOrder;
        String[] resn_threeLetterResidueCodes;
        int[] ianum_numbersOfAtomsInResidues;
        public Charges(String[][] cind_atomsByOrder,float[][] chg_chargesByOrder, String[] resn_threeLetterResidueCodes, int[] ianum_numbersOfAtomsInResidues ) {
            this.cind_atomsByOrder = cind_atomsByOrder;
            this.chg_chargesByOrder = chg_chargesByOrder;
            this.resn_threeLetterResidueCodes =resn_threeLetterResidueCodes;
            this.ianum_numbersOfAtomsInResidues = ianum_numbersOfAtomsInResidues;
        }
}
