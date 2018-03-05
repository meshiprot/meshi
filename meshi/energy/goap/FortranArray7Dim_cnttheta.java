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

        public void setAll(float value) {
            for (int i0 = 0; i0 < cnttheta.length; i0++) {
                float[][][][][][] a0 = cnttheta[i0];
                for (int i1 = 0; i1 < a0.length; i1++) {
                    float[][][][][] a1 = a0[i1];
                    for (int i2 = 0; i2 < a1.length; i2++) {
                        float[][][][] a2 = a1[i2];
                        for (int i3 = 0; i3 < a2.length; i3++) {
                            float[][][] a3 = a2[i3];
                            for (int i4 = 0; i4 < a3.length; i4++) {
                                float[][] a4 = a3[i4];
                                for (int i5 = 0; i5 < a4.length; i5++) {
                                    float[] a5 = a4[i5];
                                    for (int i6 = 0; i6 < a5.length; i6++) {
                                        a5[i6] = value;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

    public void replaceAll(float from, float to) {
        for (int i0 = 0; i0 < cnttheta.length; i0++) {
            float[][][][][][] a0 = cnttheta[i0];
            for (int i1 = 0; i1 < a0.length; i1++) {
                float[][][][][] a1 = a0[i1];
                for (int i2 = 0; i2 < a1.length; i2++) {
                    float[][][][] a2 = a1[i2];
                    for (int i3 = 0; i3 < a2.length; i3++) {
                        float[][][] a3 = a2[i3];
                        for (int i4 = 0; i4 < a3.length; i4++) {
                            float[][] a4 = a3[i4];
                            for (int i5 = 0; i5 < a4.length; i5++) {
                                float[] a5 = a4[i5];
                                for (int i6 = 0; i6 < a5.length; i6++) {
                                    if (a5[i6] == from)
                                        a5[i6] = to;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    public void trim(float value) {
        for (int i0 = 0; i0 < cnttheta.length; i0++) {
            float[][][][][][] a0 = cnttheta[i0];
            for (int i1 = 0; i1 < a0.length; i1++) {
                float[][][][][] a1 = a0[i1];
                for (int i2 = 0; i2 < a1.length; i2++) {
                    float[][][][] a2 = a1[i2];
                    for (int i3 = 0; i3 < a2.length; i3++) {
                        float[][][] a3 = a2[i3];
                        for (int i4 = 0; i4 < a3.length; i4++) {
                            float[][] a4 = a3[i4];
                            for (int i5 = 0; i5 < a4.length; i5++) {
                                boolean found = false;
                                float[] a5 = a4[i5];
                                for (int i6 = 0; i6 < a5.length; i6++) {
                                    if (a5[i6] != value)
                                        found = true;
                                }
                                if (!found)
                                    a4[i5] = null;
                            }
                        }
                    }
                }
            }
        }
    }


}
