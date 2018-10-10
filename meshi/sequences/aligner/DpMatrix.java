/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.sequences.aligner;

import java.io.IOException;

import meshi.sequences.*;
import meshi.util.Utils;

/**
 * Created by IntelliJ IDEA.
 * User: Roy
 * Date: 29/07/2005
 * Time: 15:47:31
 * To change this template use File | Settings | File Templates.
 */
public class DpMatrix {
    protected MeshiSequence sequence1, sequence2; // One dimensional sequenceAlignments
	private double minScore;
    public final CellScorer cellScorer;
    public final static int asciiA=65;
    public final static int UP=0, LEFT=1, DIAGONAL=2;
    private Cell[][] cellMatrix;
    private Cell[][][] bestRoutesMatrix;
    private double[][][] scoresMatrix;


    public DpMatrix(MeshiSequence sequence1, MeshiSequence sequence2, CellScorer scrr, double minScore) {

        this.sequence1 = sequence1;
        this.sequence2 = sequence2;
		this.minScore = minScore;
        cellScorer = scrr;
        cellMatrix = new Cell[sequence1.size() + 2][sequence2.size() + 2];
        bestRoutesMatrix= new Cell[sequence1.size()+ 2][sequence2.size()+ 2][3];
        scoresMatrix= new double[sequence1.size() + 2][sequence2.size()+ 2][5];
        for (int row = 0; row < sequence1.size() + 2; row++) {
            for (int col = 0; col < sequence2.size() + 2; col++) {
                Cell newCell = new Cell(row, col, this);
                setCell(row, col, newCell);
                /*if (matrix[row][col].maxScore<10&&row<20&&col<20)
                	System.out.print(" "+matrix[row][col].maxScore+" ");//added by Tommer 15.9.14
                else
                	if (row<20&&col<20)
                	System.out.print(matrix[row][col].maxScore+" ");*/
               
            }
            //if (row<20)
            	//System.out.print("\n");//added by Tommer 15.9.14
        }
    }

    public void setCell(int row, int column, Cell cell) {
        cellMatrix[row][column] = cell;
    }

    public Cell getCell(int rowNumber, int colNumber) {
    	if ((rowNumber < 0) | (rowNumber > sequence1.size() + 1))
    		return null;
		if ((colNumber < 0) | (colNumber > sequence2.size() + 1))
			return null;
        return cellMatrix[rowNumber][colNumber];
    }

	public SequenceAlignment backTrack(ResidueAlignmentMethod method) {
		//sequence1.printAsLine();//check, Tommer
		//sequence2.printAsLine();

		int maxI = cellMatrix.length - 1;
		int maxJ = cellMatrix[0].length - 1;
		try {
			return backTrack(cellMatrix[maxI][maxJ], minScore);
//			return backTrack(cellMatrix[cellMatrix.length-1][cellMatrix[0].length-1]);
		} catch (Exception ex) {
			ex.printStackTrace();
			throw new RuntimeException(ex);
		}

	}

	public SequenceAlignment backTrack(Cell cell, double minScore) {
        SequenceAlignment inverseAlignment;
        inverseAlignment=DpMatrix.inverseAlignment(cell, minScore, sequence1, sequence2);
 		return DpMatrix.reverseAlignment(inverseAlignment);

    }

    //OLD METHOD, replaced by Tommer 21.9.14
    /*public char rowChar(int index) {
        if (index == 0) return SequenceAlignmentCell.GAP_CHAR;
        return ((SequenceAlignmentCell) sequence1.get(index - 1).cell(0)).getChar();
    }*/
    
  //NEW METHOD, edited by Tommer 21.9.14
    public char rowChar(int index) {
        return ((SequenceAlignmentCell) sequence1.get(index - 1).cell(0)).getChar();
    }

    //OLD METHOD, replaced by Tommer 21.9.14
    /*public char columnChar(int index) {
        if (index == 0) return SequenceAlignmentCell.GAP_CHAR;
        return ((SequenceAlignmentCell) sequence2.get(index - 1).cell(0)).getChar();
    }*/
    
  //NEW METHOD, edited by Tommer 21.9.14
    public char columnChar(int index) {
        return ((SequenceAlignmentCell) sequence2.get(index - 1).cell(0)).getChar();
    }

    
	public static double backTrackScore(String sequences, AlignmentScheme alSchem) throws IOException{//this method was intended for debugging purposes, Tommer 30.9.14
		int inGap=0, inGapTemp;
		double scoreBackTrack=0;
		String chains[]=new String[2];
		sequences=sequences.substring(sequences.indexOf('\n')+1);
		chains[0]=sequences.substring(0,sequences.indexOf('\n'));
		chains[1]=sequences.substring(sequences.indexOf('\n')+1);
		for(int i=0;i<chains[0].length()&&chains[0].charAt(i)!='*';i++){
			if(chains[0].charAt(i)!='-'&&chains[1].charAt(i)!='-'
					&&chains[0].charAt(i)!=chains[1].charAt(i)){
				inGap=0;
				scoreBackTrack+=alSchem.substitutionMatrix()
						[alSchem.letterToIndex()[chains[0].charAt(i)-asciiA]][alSchem.letterToIndex()[chains[1].charAt(i)-asciiA]];//for score count
			}
			else{
				if(chains[0].charAt(i)=='-'||chains[1].charAt(i)=='-'){
					
					scoreBackTrack+=alSchem.inGapPenalty();//Tommer 30.9.14
					inGapTemp=inGap;
					if(chains[0].charAt(i)=='-'){
						inGap=1;//gap in first sequence
					}
					else{
						inGap=2;//gap in second sequence
					}
					if(inGap!=inGapTemp){
						scoreBackTrack+=alSchem.gapStartPenalty()-alSchem.inGapPenalty();//Tommer 30.9.14
					}
				}
				else {
					inGap = 0;
					scoreBackTrack += alSchem.substitutionMatrix()[alSchem.letterToIndex()[chains[0]
							.charAt(i)-asciiA]][alSchem.letterToIndex()[chains[1].charAt(i)-asciiA]];// for score count
				}
			}		
		}
		return scoreBackTrack;
		
	}
	
	public Cell[][][] bestRoutesMatrix(){
		return bestRoutesMatrix;
	}
	
	public double[][][] scoresMatrix(){
		return scoresMatrix;
	}
	
	private static SequenceAlignment inverseAlignment(Cell from, double minScore, MeshiSequence sequence1, MeshiSequence sequence2){
		Cell back = from.back;
		Cell prev = null;
		int rowNumber;
		int colNumber;
		SequenceAlignmentColumn column;
		SequenceAlignment inverseAlignment=new SequenceAlignment();
		inverseAlignment.comments.add(sequence1.comment());
        inverseAlignment.comments.add(sequence2.comment());
        inverseAlignment.setScore(from.neighborScore);

		while((from != null)) {// (from.maxScore >= minScore)){
			rowNumber = from.rowNumber;
    		colNumber = from.colNumber;
 			SequenceAlignmentCell cell1, cell2;
			if (((rowNumber == 0) | (rowNumber == sequence1.size() + 1)) ||
					((from == prev.leftCell) ))
				cell1 = new SequenceAlignmentCell();
			else cell1 = sequence1.cell(rowNumber - 1);

			if (((colNumber == 0) | (colNumber == sequence2.size() + 1)) ||
					((from == prev.upCell) ))
				cell2 = new SequenceAlignmentCell();
			else cell2 = sequence2.cell(colNumber - 1);

			column = new SequenceAlignmentColumn(cell1, cell2);


            inverseAlignment.add(column);//name modified Tommer 9.9.14

			prev = from;
			from = back;
			if (from != null) {
				back = from.back;
			}
			else {back = null;}
        }

        return inverseAlignment;

	}
	
	private static SequenceAlignment reverseAlignment(SequenceAlignment prototype){
		SequenceAlignment out= new SequenceAlignment();
		 for (int i = prototype.size() - 1; i >= 0; i--)//tempInvAlign replaces inverseAlignment Tommer 9.9.14
	            out.add(prototype.get(i));
		 out.comments.add(prototype.comments.get(0));
	     out.comments.add(prototype.comments.get(1));
	     out.setScore(prototype.score());
	     return out;
	}

	private void printMax() {
		for (int i = 0; i < cellMatrix.length; i++){
			for (int j = 0; j < cellMatrix[0].length; j++) {
				System.out.print(cellMatrix[i][j].neighborScore + cellMatrix[i][j].internalScore + " ");
			}
			System.out.println();
		}
	}
	
}


   

	
