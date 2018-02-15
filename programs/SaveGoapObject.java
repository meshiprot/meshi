package programs;

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
    }
    public static void save() throws Exception {

        GoapCreator goapCreator = new GoapCreator();
        goapCreator.parametersPath="D:\\workspace\\meshi\\parameters\\meshiPotential\\GOAP\\";
        Fort31 fort31 = goapCreator.getFort31();
        System.out.println(fort31);
        FileOutputStream fos = new
                FileOutputStream("fort31_db");
        GZIPOutputStream gz = new GZIPOutputStream(fos);
        ObjectOutputStream oos = new ObjectOutputStream(gz);
        oos.writeObject(fort31);
        oos.flush();
        oos.close();
        fos.close();
    }
    public static void load() throws Exception{

        FileInputStream fis = new FileInputStream("fort31_db");
        GZIPInputStream gs = new GZIPInputStream(fis);
        ObjectInputStream ois = new ObjectInputStream(gs);
        Fort31 fort31 = (Fort31) ois.readObject();
        System.out.println(fort31);

        ois.close();
        fis.close();
    }
}
