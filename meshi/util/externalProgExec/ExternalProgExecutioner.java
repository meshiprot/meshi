package meshi.util.externalProgExec;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Created by user on 08/03/2018.
 */
public abstract class ExternalProgExecutioner {

    public ExternalProgExecutioner(){    }

    public abstract String exec(String[] args) throws IOException, InterruptedException;

    protected String execProg(String[] executionString)throws IOException,InterruptedException{
        ProcessBuilder pb = new ProcessBuilder(executionString);
        Process process = pb.start();
        int errCode = process.waitFor();
        return output(process.getInputStream());
    }
    protected String output(InputStream inputStream) throws IOException {
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
}
