package programs;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Created by user on 05/03/2018.
 */
public class executeDssp {
    public static void main(String[] args) throws InterruptedException,
            IOException {
        ProcessBuilder pb = new ProcessBuilder(args);
        System.out.println("Run command");
        Process process = pb.start();
        int errCode = process.waitFor();
        System.out.println("Command executed, any errors? " + (errCode == 0 ? "No" : "Yes"));
        System.out.println("Output:\n" + output(process.getInputStream()));
    }

    private static String output(InputStream inputStream) throws IOException {
        StringBuilder sb = new StringBuilder();
        BufferedReader br = null;
        try {
            br = new BufferedReader(new InputStreamReader(inputStream));
            String line = null;
            while ((line = br.readLine()) != null) {
                sb.append(line + System.getProperty("line.separator"));
            }
        } finally {
            br.close();
        }
        return sb.toString();
    }
}
