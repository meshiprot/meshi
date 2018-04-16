package meshi.applications.contactMap;

import meshi.molecularElements.Residue;
import meshi.parameters.ResidueType;

import java.io.Serializable;

public class ContactMapResidue implements Serializable{
    public final ResidueType type;
    public final String proteinName;
    public final int number;
    public final String chain;

    public ContactMapResidue(String proteinName, Residue residue) {
        type = residue.type();
        this.proteinName = proteinName;
        number = residue.number();
        chain = residue.chain();
    }

    public boolean equals(Object obj) {
        if (!(obj instanceof ContactMapResidue)) return false;
        ContactMapResidue other = (ContactMapResidue) obj;
        if (((proteinName != null) & (other.proteinName != null)) && (!proteinName.equals(other.proteinName)))
            return false;
        if (number != other.number) return false;
        if (!chain.equals(other.chain)) return false;
        if (type != other.type) return false;
        return true;
    }
}
