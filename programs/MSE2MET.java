package programs;

import meshi.util.file.MeshiLineReader;
import meshi.util.file.MeshiWriter;

import java.io.IOException;

/**
 * Created by chen on 02/03/2017.
 */
public class MSE2MET {
    public static void main(String[] args) throws IOException {
        MeshiLineReader reader = new MeshiLineReader(args[0]);
        MeshiWriter writer = new MeshiWriter(args[1]);
        String line = "";
        while ((line = reader.readLine()) != null) {
            if (line.startsWith("ATOM")) writer.println(line);
            else if (line.startsWith("HETATM")) {
                int index = line.indexOf("MSE");
                if (index > 0) {
                    String newLine = "ATOM  "+line.substring(6,index) + "MET" + line.substring(index + 3);
                    index = line.indexOf(" SE ");
                    if (index < 0) writer.println(newLine);
                    else {
                        newLine = "ATOM  "+line.substring(6,index) + "  SD  MET" + line.substring(index + 9);
                        writer.println(newLine);
                    }
                }
            }
        }
        writer.close();
    }
}
