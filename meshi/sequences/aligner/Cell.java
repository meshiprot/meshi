/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 29/07/2005
 * Time: 15:39:55
 * To change this template use File | Settings | File Templates.
 */
public class Cell {
    protected enum Type {NORMAL, TOP, BOTTOM, LEFT, RIGHT, FIRST, LAST}
    public final static int UP=0, LEFT=1, DIAGONAL=2, MAX=3, INTERNAL = 4;
    public final int rowNumber;
    public final int colNumber;
    public final Cell upCell, diagonalCell, leftCell;//the next cell will be either the cell in the right of this cell or the cell in bottom of this one
    public final DpMatrix dpMatrix;
    public final Type type;

    protected double startPenalty, inPenalty;
    protected Cell back, nextBack;//nextBack added by Tommer 4.9.14
    protected double neighborScore, internalScore;

    public Cell(int rowNumber, int colNumber, DpMatrix mat) {
        this.rowNumber = rowNumber;
        this.colNumber = colNumber;
        dpMatrix = mat;
        type = getType(rowNumber, colNumber);
        upCell       = dpMatrix.getCell(rowNumber - 1, colNumber);
        leftCell     = dpMatrix.getCell(rowNumber, colNumber - 1);
        diagonalCell = dpMatrix.getCell(rowNumber - 1, colNumber -1);
        //the score for the cell
        dpMatrix.cellScorer.getScores(this);
    }

    private Type getType(int row, int column) {
        if (row == 0) {
            if (column == 0) return Type.FIRST;
            else return Type.TOP;
        }
        if (column == 0) return Type.LEFT;
        if (row == dpMatrix.sequence1.size() + 1) {
            if (column == dpMatrix.sequence2.size() + 1)
                return Type.LAST;
            else return Type.BOTTOM;
        }
        if (column == dpMatrix.sequence2.size() + 1)
            return Type.RIGHT;
        return Type.NORMAL;
    }

    public void setBack(Cell back, Cell nextBack) {//modified Tommer 4.9.14
//        if ((back != null) &&(nextBack != back.back))
//            throw new RuntimeException("This is weird "+back+" ; "+back.back+" ; " + nextBack);
        this.back = back;
        this.nextBack=nextBack;
    }
    

    public Cell getBack() {
        return back;
    }

    public String toString() {
        String backString = "";
        if (back != null)
            backString += " ("+back.rowNumber+" , "+back.colNumber+")";
        else backString = " (null)";
        return "cell " + rowNumber + " " + colNumber + " " +
                neighborScore + " " + internalScore + " " + (neighborScore + internalScore) + backString;
    }

}
