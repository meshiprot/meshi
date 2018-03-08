package meshi.util.externalProgExec;

import meshi.util.CommandList;

import java.io.IOException;

/**
 * Created by user on 08/03/2018.
 */
public class ScwrlExec extends ExternalProgExecutioner {

    private static final String commandsKeyName = "scwrlExecPath";
    String[] executionString;

    String progPath;

    public ScwrlExec(CommandList commands){
        super();
        progPath=commands.firstWord(commandsKeyName).secondWord();
    }

    /**
     * @args[0] - String pdbInputFilePath
     * @args[1] - String pdbOutputFilePath
     */
    public String exec(String args[]) throws IOException, InterruptedException {
        executionString = new String[]{progPath,"-i",args[0],"-o",args[1]};
        return execProg(executionString);
    }
}
