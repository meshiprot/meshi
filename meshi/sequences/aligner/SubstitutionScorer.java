/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

import meshi.sequences.*;
import meshi.util.Utils;
import org.omg.PortableInterceptor.INACTIVE;

public class SubstitutionScorer implements CellScorer {
	public final static int asciiA=65, numAminoAcids=20, numWildChars=3, numEnglishLetters=26;
    private double gapStartPenalty;//added/modified by Tommer, 2.9.14
    private double inGapPenalty;//modified by Tommer, 2.9.14
    private int substitutionMatrix[][]=new int[numAminoAcids+numWildChars][numAminoAcids+numWildChars];//added 1.9.14
    private int letterToIndex[] = new int[numEnglishLetters];//might need deleting 1.9.14

    
    public SubstitutionScorer(SubstitutionMatrix alSchem) {//signature changed Tommer 1.9.14
    	this.gapStartPenalty =alSchem.gapStartPenalty();
        this.inGapPenalty = alSchem.inGapPenalty();
        this.substitutionMatrix=alSchem.substitutionMatrix();//Added
        this.letterToIndex=alSchem.letterToIndex();
    }

    
    public void getScores(Cell cell) {
        DpMatrix matrix = cell.dpMatrix;

		cell.startPenalty = gapStartPenalty;
        cell.inPenalty    = inGapPenalty;
		cell.internalScore =internalScore(matrix, cell, substitutionMatrix, letterToIndex);
        if (cell.type != Cell.Type.NORMAL) {
			cell.startPenalty = cell.inPenalty = 0;
			cell.internalScore = 0;
		}



		double scoreLeft     = getScoreLeft(cell);
		double scoreUp       = getScoreUp(cell);
		double scoreDiagonal = getScoreDiagonal(cell);
		cell.neighborScore = maxMultiple(scoreLeft, scoreUp, scoreDiagonal);
		cell.neighborScore = maxMultiple(scoreLeft, scoreUp, scoreDiagonal);
		if (cell.type == Cell.Type.LEFT)
			cell.back = cell.upCell;
		else if (cell.type == Cell.Type.TOP)
			cell.back = cell.leftCell;
		else if (cell.neighborScore == scoreDiagonal) {
			cell.back = cell.diagonalCell;
		} else if (cell.neighborScore == scoreLeft) {
			cell.back = cell.leftCell;
		}else if (cell.neighborScore == scoreUp) {
			cell.back = cell.upCell;
		} else throw  new RuntimeException("This is weird");

		if (cell.back != null)
			cell.nextBack = cell.back.back;
		else cell.nextBack = null;


    }

    private double getScoreLeft(Cell cell) {
    	if ((cell.type == Cell.Type.FIRST ) | (cell.type == Cell.Type.LEFT ))
    		return 0;
    	Cell left = cell.leftCell; // cannot be null now
		double out = left.neighborScore;
		if (left.back != null) {
			if (left.back == left.leftCell)
				out += cell.inPenalty;
			else if (left.back == left.upCell)
				out += cell.startPenalty;
			else if (left.back == left.diagonalCell)
				out += cell.startPenalty;
			else throw new RuntimeException("This is weird");
		}
    	return out;
	}

	private double getScoreUp(Cell cell) {
		if ((cell.type == Cell.Type.FIRST ) | (cell.type == Cell.Type.TOP ))
			return 0;
		Cell up = cell.upCell; // cannot be null now
		double out = up.neighborScore;
		if (up.back != null) {
			if (up.back == up.leftCell)
				out += cell.startPenalty;
			else if (up.back == up.upCell)
				out += cell.inPenalty;
			else if (up.back == up.diagonalCell)
				out += cell.startPenalty;
			else throw new RuntimeException("This is weird");
		}
		return out;
	}
	private double getScoreDiagonal(Cell cell) {
		if ((cell.type == Cell.Type.FIRST ) | (cell.type == Cell.Type.TOP ) | (cell.type == Cell.Type.LEFT ))
			return 0;
		Cell diagonal = cell.diagonalCell; // cannot be null now
		double out = diagonal.neighborScore + diagonal.internalScore;
		if (diagonal.back != diagonal.diagonalCell)
			out += cell.startPenalty;
		return out;
	}

	public double maxMultiple(double d1, double d2, double d3){
		double[] arr = {d1, d2, d3};
		return maxMultiple(arr);
	}
	public double maxMultiple(double[] array){
    	try{
    	if (array==null||array.length<1){
    		throw new RuntimeException("null array or empty array"); 
    	}
    	}
    	catch(RuntimeException e){
    		System.out.println("null array or empty array");
    	}
    	
    	double temp=array[0];
    	for (int i=0;i<array.length;i++){
    		temp=Math.max(temp, array[i]);
    	}
    	return temp;
    }
    
    public static double internalScore(DpMatrix matrix, Cell cell, int[][] subMat, int[]revOrder){
    	if ((cell.rowNumber == 0) | (cell.colNumber == 0) | (cell.rowNumber == matrix.sequence1.size()+1) | (cell.colNumber == matrix.sequence2.size() + 1))
			return 0;
    	char rowChar = matrix.rowChar(cell.rowNumber);
        char columnChar = matrix.columnChar(cell.colNumber);
        
        double internalScore;
        
       
        if ((rowChar == SequenceAlignmentCell.GAP_CHAR) |
                (columnChar == SequenceAlignmentCell.GAP_CHAR) |
                (rowChar == SequenceAlignmentCell.WILDCARD_CHAR) |
                (columnChar == SequenceAlignmentCell.WILDCARD_CHAR))
        	internalScore = 0;
        //else internalScore = ((rowChar == columnChar) ? 5 :0);//Edited by Tommer 1.9.14, mismatch penalty was 0
        
        /*the following implementation replaces the above line. it should do the trick
         (of allowing for the use of a 2D matrix score) ASSUMING this is really the only
          relevant place were score is given for matches/mismatches. this is to be checked.  */
        
       else{
        	if(revOrder[(int)rowChar-asciiA]==-1){
        		System.out.println("error: unrecognized residue letter");
        		internalScore =0;
        	}
        	else
        		internalScore =subMat[revOrder[(int)rowChar-65]][revOrder[(int)columnChar-65]];
        		
       } 
        return internalScore;
    }
}

