package programsTests;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

/**
 * Created by user on 15/02/2018.
 */
public class CompressedReaderTest {
    final static int BUFFER = 2048;
    static String fileName = "D:\\workspace\\meshi\\parameters\\meshiPotential\\GOAP\\fort.zip";
    public static void main(String[] args) {
        try {
            /*InputStream fin = new FileInputStream("/data/file.bz2");
            BufferedInputStream bis = new BufferedInputStream(fin);
            CompressorInputStream input = new CompressorStreamFactory().createCompressorInputStream(bis);
            BufferedReader br = new BufferedReader(new InputStreamReader(input,"UTF-8"));

            String line = "";
            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }
            */

            BufferedOutputStream dest = null;
            FileInputStream fis = new
                    FileInputStream(fileName);
            ZipInputStream zis = new
                    ZipInputStream(new BufferedInputStream(fis));
            ZipEntry entry;
            while((entry = zis.getNextEntry()) != null) {
                System.out.println("Extracting: " +entry);
                int count;
                byte data[] = new byte[BUFFER];
                String str = "";
                while ((count = zis.read(data, 0, BUFFER)) != -1) {
                    str = new String(data);

                    System.out.println(str);
                    //dest.write(data, 0, count);
                }

                // write the files to the disk
                /*
                FileOutputStream fos = new
                        FileOutputStream(entry.getName());
                dest = new
                        BufferedOutputStream(fos, BUFFER);
                while ((count = zis.read(data, 0, BUFFER))
                        != -1) {
                    dest.write(data, 0, count);
                }
                dest.flush();
                dest.close();
                */
            }
            zis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
