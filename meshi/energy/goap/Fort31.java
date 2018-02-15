package meshi.energy.goap;

import java.io.Serializable;

/**
 * Created by chen on 05/07/2015.
 */
public class Fort31 implements Serializable{
        public FortranArray7Dim_cnttheta cnttheta_unknown;
        public int[] map_unknown;
    int ig_s_parameter;
        public Fort31(FortranArray7Dim_cnttheta cnttheta_unknown, int[] map_unknown, int ig_s_parameter) {
            this.cnttheta_unknown = cnttheta_unknown;
            this.map_unknown = map_unknown;
            this.ig_s_parameter = ig_s_parameter;
        }
}
