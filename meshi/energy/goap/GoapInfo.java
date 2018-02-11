package meshi.energy.goap;

import meshi.energy.EnergyInfoElement;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfo;

/**
 * Created by chen on 05/07/2015.
 */
public class GoapInfo extends EnergyInfoElement{
    public final MeshiInfo dFIRE, goapAG;
    public GoapInfo() {
        super(InfoType.GOAP,"GOAP");
        getChildren().add(dFIRE  = new MeshiInfo(InfoType.D_FIRE, Double.valueOf(-9999), "Orientation INdependent component of GOAP"));
        getChildren().add(goapAG = new MeshiInfo(InfoType.GOAP_AG, Double.valueOf(-9999), "Orientation  dependent component of GOAP"));
    }

}
