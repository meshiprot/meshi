package meshi.energy.goap;

import java.io.Serializable;

/**
 * Created by chen on 05/07/2015.
 */
public class FortranArray5Dim_pot implements Serializable{

        float[][][][][] potential_pot = new float[GoapCreator.ibin_unknown][20][15][20][15];//last index made first! 31.3.15

        protected float get_potential(int i1,int i2, int i3, int i4, int i5){
            return potential_pot[i5-1][i1-1][i2-1][i3-1][i4-1];//last index made first! 31.3.15
        }
        protected void set_potential(int i1,int i2, int i3, int i4, int i5, float value){
            potential_pot[i5-1][i1-1][i2-1][i3-1][i4-1]=value;//last index made first! 31.3.15
        }
}
