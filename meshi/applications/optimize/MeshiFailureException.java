package meshi.applications.optimize;

import meshi.util.file.MeshiWriter;

public class MeshiFailureException extends RuntimeException{
    public MeshiFailureException(String message) {
        super(message);
    }
}
