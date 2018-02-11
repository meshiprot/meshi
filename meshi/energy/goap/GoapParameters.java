package meshi.energy.goap;

import meshi.energy.Parameters;

/**
 * Created by chen on 05/07/2015.
 */
public class GoapParameters implements Parameters{
    protected Fort21 fort21;
    protected Fort31 fort31;
    protected Charges charges;
    protected SideGeometry sideGeometry;

    public GoapParameters(Fort21 fort21, Fort31 fort31, Charges charges, SideGeometry sideGeometry) {
        this.fort21 = fort21;
        this.fort31 = fort31;
        this.charges = charges;
        this.sideGeometry=sideGeometry;
    }
}
