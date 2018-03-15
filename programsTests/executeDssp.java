package programsTests;

import java.io.*;

/**
 *  Created by user on 05/03/2018.
 *
 *   */
public class executeDssp {
    public static void main(String[] args) throws InterruptedException,
            IOException {
        ProcessBuilder pb = new ProcessBuilder(args);

        pb.redirectErrorStream(true);
        Process process = pb.start();
        OutputStream pos = process.getOutputStream();

        readModels("/home/cluster/users/siditom/data/Amarda/amarda_test/2HG6_rosettaDecoys.test.pdb",0,1,pos);
        pos.close();

        int errCode = process.waitFor();
        System.out.println("Output:\n" + output(process.getInputStream()));
    }

    private static String output(InputStream inputStream) throws IOException {
        StringBuilder sb = new StringBuilder();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new InputStreamReader(inputStream));
            String line = null;

            while ((line = br.readLine()) != null) {
                //	System.out.println(line);
                sb.append(line + System.getProperty("line.separator"));
            }
        }catch (Exception e){
            System.out.println(e);
        } finally {
            br.close();
        }
        return sb.toString();
    }
    private static void readModels(String filepath,long inxStartModel,int numOfModels,OutputStream pos) {
        try {

            RandomAccessFile raf = new RandomAccessFile(filepath, "r");//Open our file with read/write access
            StringBuilder sb = new StringBuilder();

            raf.seek(inxStartModel);
            String line = "";
            int iModel = 0;
            byte[] b=null;
            while (line != null && iModel < numOfModels) {
                line = raf.readLine();
                while (line != null && !line.contains("END MODEL")) {
                    sb.append(line + System.getProperty("line.separator"));
                    b = line.getBytes();
                    line = raf.readLine();
                }
                if (line !=null) {
                    sb.append(line+ System.getProperty("line.separator"));
                }
                iModel++;
            }
            System.out.println(sb.toString());
            pos.write(sb.toString().getBytes());
            raf.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
