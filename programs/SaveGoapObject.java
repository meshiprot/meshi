package programs;

import meshi.energy.goap.Fort21;
import meshi.energy.goap.Fort31;
import meshi.energy.goap.GoapCreator;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Created by user on 15/02/2018.
 */
public class SaveGoapObject {
    public static void main(String[] argv) throws Exception {
        SaveGoapObject.save();
        SaveGoapObject.load();
        //float[][][][][][][] cnttheta=new float[20][15][20][15][GoapCreator.ibin_unknown][5][12];
        //System.out.println("Hello");
        //System.out.println("Hello");

    }
    public static void save() throws Exception {

        GoapCreator goapCreator = new GoapCreator();
        //goapCreator.parametersPath="D:\\workspace\\meshi\\parameters\\meshiPotential\\GOAP\\";
        goapCreator.parametersPath="D:\\workspace\\meshi\\parameters_backup\\";
        Fort21 fort21 = goapCreator.getFort21();
        System.out.println(fort21);
        FileOutputStream fos = new
                FileOutputStream("parameters/meshiPotential/GOAP/fort21_db");
        GZIPOutputStream gz = new GZIPOutputStream(fos);
        ObjectOutputStream oos = new ObjectOutputStream(gz);
        oos.writeObject(fort21);
        oos.flush();
        oos.close();
        fos.close();
    }
    public static void load() throws Exception{

        FileInputStream fis = new FileInputStream("D:\\workspace\\meshi\\parameters\\meshiPotential\\GOAP\\fort31_db");
        GZIPInputStream gs = new GZIPInputStream(fis);
        ObjectInputStream ois = new ObjectInputStream(gs);
        Fort31 fort31 = (Fort31) ois.readObject();
        System.out.println(fort31);

        ois.close();
        fis.close();
    }
}
