package meshi.util.externalProgExec;

import meshi.util.CommandList;

import java.io.IOException;

/**
 * Created by user on 08/03/2018.
 */
public class DsspExec extends ExternalProgExecutioner {

    private static final String commandsKeyName = "dsspExecPath";
    String[] executionString;

    String progPath;

    public DsspExec(CommandList commands){
        super();
        progPath=commands.firstWord(commandsKeyName).secondWord();
    }

    /**
     * @args[0] - String pdbInputFilePath
     * @args[1] - String dsspOutputFilePath
     */
    public String exec(String args[]) throws IOException, InterruptedException {
        executionString = new String[]{progPath,args[0],args[1]};
        return execProg(executionString);
    }
}
