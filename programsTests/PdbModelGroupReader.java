package programsTests;

import meshi.util.file.MeshiWriter;

import java.io.*;
import java.util.Vector;

/**
 * Created by user on 28/02/2018.
 */
public class PdbModelGroupReader {


    public static Vector<Long> modelstart = new Vector<Long>();
    public static int numOfModelsInGroup;
    public static int totalNumOfModels;
    public static String filepath;

    public static void main(String[] args) throws Exception{
        filepath="D:\\workspace\\meshi\\siditom_tests\\amarda_test\\2HG6_rosettaDecoys.test.pdb";
        totalNumOfModels = 0;
        numOfModelsInGroup = 1;
        getModelInxs(filepath,numOfModelsInGroup);
        writeInxs();
        //for (long l : modelstart) System.out.println(l);
        //readModels(filepath,3,2);
    }

    public static void writeInxs() throws IOException {
        MeshiWriter mw = new MeshiWriter(new File(filepath+".inx"));
        mw.println(filepath);
        mw.print(totalNumOfModels);
        int imodel=1;
        for (long inx : modelstart){
            mw.print("\n"+inx+" "+imodel+" "+numOfModelsInGroup);
            imodel += numOfModelsInGroup;
        }
        mw.close();
    }

    /**
     * Deprecated - example how to read a file using RandomAccessFile.
     *
     * @param filepath
     * @param inxStartModel
     * @param numOfModels
     */
    private static void readModels(String filepath,int inxStartModel,int numOfModels) {
        try {

            RandomAccessFile raf = new RandomAccessFile(filepath, "r");//Open our file with read/write access

            raf.seek(modelstart.elementAt(inxStartModel-1).longValue());
            String line = "";
            int iModel = 0;
            while (line != null && iModel < numOfModels) {
                line = raf.readLine();
                while (line != null && !line.contains("END MODEL")) {
                    System.out.println(line);
                    line = raf.readLine();
                }
                if (line !=null) System.out.println(line);
                iModel++;
            }

            raf.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    public static void getModelInxs(String filepath, int numOfModels) {
        try {
            FileInputStream fis = new FileInputStream(filepath);

            DataInputStream dis = new DataInputStream(fis);
            BufferedInputStream bis = new BufferedInputStream(dis);
            modelstart.add((long)0);
            int iModel = 0;
            int nb = 0;
            int loc = -1;
            byte[] b = new byte[512];
            long nbytes = 0;
            String str = "";
            String remainerStr = "";
            while(dis.available()>0) {
                nb = dis.read(b);
                str = remainerStr + new String(b,0,nb);
                remainerStr = str.substring(str.lastIndexOf("\n"));
                str = str.substring(0,str.lastIndexOf("\n"));
                loc = str.indexOf("END MODEL");
                if (loc > 0){
                    iModel++;
                    totalNumOfModels++;
                    if (iModel == numOfModels) {
                        modelstart.add(nbytes + str.substring(0, loc).getBytes().length + "END MODEL\n".getBytes().length);
                        iModel = 0;
                    }
                }
                nbytes += str.getBytes().length;
            }
            dis.close();
            fis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    public static void getModelEnds(String filepath){
        try {
            RandomAccessFile raf = new RandomAccessFile(filepath, "r");//Open our file with read/write access
            String line = "";
            line = raf.readLine();
            while (line != null) {
                if (line.contains("END MODEL")){
                    modelstart.add(raf.getFilePointer());
                }
                line = raf.readLine();
            }
            raf.close();//Close our filestream.
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
