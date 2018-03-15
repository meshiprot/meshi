package programs;

import meshi.util.MeshiProgram;
import meshi.util.file.MeshiWriter;

import java.io.*;
import java.util.Vector;

/**
 * Created by user on 28/02/2018.
 */
public class InitPdbIndex extends MeshiProgram{

    public static Vector<Long> modelstart = new Vector<Long>();
    public static int totalNumOfModels;
    public static String inFileName;
    public static int nModels;
    public static String outPath;

    public static void main(String[] args) throws Exception{
        init(args);
        totalNumOfModels = 0;
        getModelInxs(inFileName,1);
        writeInxs();
    }

    public static long[] getGroup(String inxFile, int iGroup) throws IOException{
            RandomAccessFile raf = new RandomAccessFile(inxFile, "r");
            String line = raf.readLine();
            while (line != null && line.compareTo("--INDEX--")!=0) line=raf.readLine();
            for (int i=0; line != null && i<iGroup; i++){
                line = raf.readLine();
            }
            if (line == null) throw new RuntimeException("Error - incorrect group index - iGroup > Number of Groups in PDB file.");

            String[] indicesSTR = line.split(",");


            long[] indices = new long[indicesSTR.length];
            for (int i=0; i < indices.length; i++){
                indices[i] = Long.parseLong(indicesSTR[i]);
            }
            return indices;
    }

    public static void writeInxs() throws IOException {
        MeshiWriter mw = new MeshiWriter(new File(inFileName+".inx"));
        mw.println("--HEADER--");
        mw.println("PDB file path: "+inFileName);
        mw.println("Total Number Of Models: "+totalNumOfModels);
        mw.println("Number Of Models In A Group: "+nModels);
        mw.println("Number Of Groups In PDB File: "+Math.ceil(((double)totalNumOfModels)/nModels));
        mw.println("\n--INDEX--");
        for (int i=0; i < modelstart.size()-1; i++){

            if ((i+1)%nModels == 0 || i == totalNumOfModels-1) mw.println(modelstart.get(i));
            else mw.print(modelstart.get(i)+",");
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

    private static void init(String[] argv) {
        //                 0            1            2            3                4              5
        String[] keys = {"inFileName","nModels","out"};
        String[] arguments = getArguments(keys, argv);

        inFileName = arguments[0];
        nModels = Integer.parseInt(arguments[1]);
        outPath=arguments[2];

    }
}
