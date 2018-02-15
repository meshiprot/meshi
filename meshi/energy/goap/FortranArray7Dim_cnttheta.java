package meshi.energy.goap;

import java.io.Serializable;

/**
 * Created by chen on 05/07/2015.
 */
public class FortranArray7Dim_cnttheta implements Serializable{

        float[][][][][][][] cnttheta=new float[20][15][20][15][GoapCreator.ibin_unknown][5][12];//last 2 indices switched! 31.3.15
        protected float
        get_cnttheta(int i1,int i2, int i3, int i4, int i5, int i6, int i7){
            return cnttheta[i1-1][i2-1][i3-1][i4-1][i5-1][i7-1][i6-1];//last 2 indices switched! 31.3.15
        }
        protected void
        set_cnttheta(int i1,int i2, int i3, int i4, int i5, int i6, int i7, float value){
            cnttheta[i1-1][i2-1][i3-1][i4-1][i5-1][i7-1][i6-1]=value;//last 2 indices switched! 31.3.15
        }
}
