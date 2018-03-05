package meshi.util.externalProgExec;

import meshi.util.CommandList;

import java.io.*;

/**
 * Created by user on 05/03/2018.
 */
public class ExternalFeatureExtractor {
    private static final String[] execs = {"dsspExecPath"};
    private static final String[] args = {" --"};


    public static String[] getExternals(CommandList commands, String pdbFile) throws InterruptedException, IOException {
        String[] outputs = new String[execs.length];
        for (int i = 0; i < execs.length; i++) {
            String executionString = commands.firstWord(execs[i]).secondWord() + " " + args[i];
            String res = execProg(executionString,pdbFile);
        }
        return outputs;
    }
    public static String execProg(String executionString, String pdbFile)throws IOException,InterruptedException{
        ProcessBuilder pb = new ProcessBuilder(executionString);
        pb.redirectErrorStream(true);
        Process process = pb.start();
        OutputStream pos = process.getOutputStream();
        pos.write(pdbFile.toString().getBytes());
        pos.close();
        int errCode = process.waitFor();
        return output(process.getInputStream());
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
}
