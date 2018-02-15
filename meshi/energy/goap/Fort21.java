package meshi.energy.goap;

import java.io.Serializable;

/**
 * Created by chen on 05/07/2015.
 */
public class Fort21 implements Serializable{
        public int[] map_unknown;
        public FortranArray5Dim_pot pot_potential;
        public Fort21(int[] map_unknown, FortranArray5Dim_pot pot_potential) {
            this.map_unknown   = map_unknown;
            this.pot_potential = pot_potential;
        }
}
