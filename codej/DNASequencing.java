import java.util.Arrays;
import java.util.Map;
import java.util.ArrayList;
import java.util.TreeMap;

/**
 * DNASequencing class - all methods for sequencing
 * 
 * @author Shiv
 * @version 1.0
 */

public class DNASequencing {
	/**
	 * dbg levels -1,0,1 (comparator is dbg >. so dbg = -1 shuts all dbg messages) 
	 */
	int dbg = -1;
	
	/**
	 * chrLSeq contains the total number of nucleotides in the chromatid specified by its index. 
	 * chrLnSeq contains the total number of bytes for this chromatid
	 *          (both indexes are 0 based, so chr 1 value will be in index 0 and so on)
	 *          Chromatids 25 and 26 are toy chromatids used for testing
	 */
	private static final int[] chrLSeq = {248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,133797422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,156040895,57227415,18,18};
	private static final int[] chrLnSeq= {62239106,60548383,49573890,47553639,45384565,42701495,39836494,36284659,34598680,33449356,33771656,33318828,28591082,26760930,25497798,22584587,20814361,20093322,14654404,16111042,11677496,12704617,39010224,14306854,5,5};	
	
	private static final byte[] g_vals = {0b00000000,0b00000001,0b00000010,0b00000011};
	
	/**
	 * Constants for the chromatids used in the diferent test difficulty levels
	 */
	int[][] TDChr = {{20},
				{1,11,20},
				{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24},
				{25},
				{25,26}};
	int testDifficulty = 0;
	
	/**
	 * chrIdA   Array of the chromatid numbers used in the diferent test difficulty levels
	 * chrIdPtr keeps track of the current position where the next chromatidID should be set
	 *          Used: Both of these are used when the chromatids are passed in through passReferenceGenome and then used as a level of indirection
	 */	
	int[] chrIdA; //chr nos that are used
	int chrIdPtr = 0;
	

	/**
	 * chrNA   Type: Array with length
	 *         Array of the NN blocks collected for all chromatids
	 *         The length of NN Blocks is stored as the 0th element 
	 *         If there are 4 NN Blocks then this will be 8 and total length of the array will be 9
	 *         Each NN Block is stored as a pair in two contiguous elements, the 1st denoting start location and 2nd denoting end location
  	 *         These locations are the actual number and not the byte (B) offsets
	 * chrnCtr keeps track of the current position where the next NN block should be set
	 *         Used: convChrToB fills them up and then they are used for hashing
	 *               If the NN Blocks are stored in a file (preprocessed) then they are set by the Set functions
	 */		
	int[] chrNA; //chr N array for each chr in the order of chrIDA
	int chrNCtr = 1;

	/**
	 * chrAA   chromatid sequence stored in B format (4 chromatids to a byte)
	 *         All Ns are stored as A with chrNA keeping track of these locations
	 *         Used: convChrToB fills it up
	 *               If the NN Blocks are stored in a file (preprocessed) then they are set by the Set functions
	 */			
	byte[] chrAA ; 

	/**
	 * chrIdOffset   chromatid start lengths - Gap between chromosomes is specified by GAPBETCHR
 	 *               Used: For the output this is used to calculate the read offset in the chromatid
	 */				
	int[] chrIdOffset; 

	/**
	 * chrNAChr   Type: Array with length
	 *            An offset into the chrNA array for all chromatids
	 *            There are a total of 2 * number of chromatids used + 1 (for array length) elements in this Array
	 *            Each (startpos for chromatid,length of NN Blocks) is stored as a pair in two contiguous elements
	 */
	int[] chrNAChr; //chr N array for each chr in the order of chrIDA
	
	/**
	 * totLength total length of all the loaded chromatids in bytes (4 chr to a byte)
	 */	
	int totLength; //tot lenght of the chr in bytes (4 chr to a byte)

	/**
	 * Used for getAlignment and calculating probables
	 */
	int[] numProbablesM1 = new int[splitA.length * 2];
	int[] numProbablesM2 = new int[splitA.length * 2];
	Probable[][] probablesM1 = new Probable[splitA.length * 2][NUMPROBABLES];
	Probable[][] probablesM2 = new Probable[splitA.length * 2][NUMPROBABLES];
	
	//Hashing
	HashIndex hInd = new HashIndex();
	
	/**
	 * CONSTANTS
	 */
	static final int GAPBETCHRBYTE = 250; //The Gap between chromosomes in compressed format (4 chr to a byte)- ie 1000 numbers
	static final int GAPBETCHR = 1000; //The gap inserted between chromosomes in ChrAA
    static final int k = 45;  //Length of the string used for hashing
	
	static final int IDEALREADGAP = 450;  //Ideally the reads will be 450 apart
	static final int READLENGTH = 150;    //input read lengths for alignment

	static final int MAXREADGAP = 855; //650; //used for alignment
	static final int MINREADGAP = 278; //250; //used for alignment
	static final int MINFREQ = 3;
	
	static final int READGAPDEN = (MAXREADGAP - MINREADGAP) / 2;  //Read Gap Denominator
	
	static final int NUMPROBABLES = 10;  //number of probables to consider for alignment
	
	static final int[] splitA = {0,15,30,45,60,75,90,105}; //start positions for the string split used for hashing
	
	//For TESTING
	//int k = 4;
	//int[] splitA = {0,2};
	//static final int READLENGTH = 6;	

	/** 
	 * print function for chrNA
	 * chrNA is a single dim Array with all NN Blocks
	 * @param start  - the index at which to start printing
	 * @param select - if select = 0 , print all NN blocks, else print upto select entry
	 *
	 * @return none
	 */
	void printChrNA(int start,int select) {	
		System.out.println("chr NN Blocks length=" + chrNA[0] + " chrNCtr = " + chrNCtr);
		if (select > 0) {
			for (int j=start; j < start + select; j+=2)				
				System.out.println("start=" + Integer.toUnsignedString(chrNA[j]) + " End=" + Integer.toUnsignedString(chrNA[j+1]) + " diff="+ (chrNA[j+1] - chrNA[j]));
		}else {
			System.out.println("Differences");
			for (int j=1; j < chrNA[0] + 1; j+=2)
				if ((chrNA[j+1] - chrNA[j]) > 0)
					System.out.println("start=" + Integer.toUnsignedString(chrNA[j]) + " End=" + Integer.toUnsignedString(chrNA[j+1]) + " diff="+ (chrNA[j+1] - chrNA[j]));
		}
	}

	/** 
	 * print function for chrNAChr
	 * @param none
	 *
	 * @return none
	 */	
	void printChrNAChr() {	
		System.out.println("chrNAChr Blocks length=" + chrNAChr[0]);
		for (int j=1; j < chrNAChr[0]; j+=2)
			System.out.println("offset= " + Integer.toUnsignedString(chrNAChr[j]) + " Length= " + Integer.toUnsignedString(chrNAChr[j+1]));
	}

	/** 
	 * print function for chrIDOffset
	 * @param testDifficulty 
	 *
	 * @return none
	 */		
	void printChrIdOffset(int testDifficulty) {	
		System.out.println("**ChrIdOffset**" );
		for (int x : chrIdOffset)
			System.out.print(Integer.toUnsignedString(x) + ",");
		System.out.println();
		
		int lSeq = chrLSeq[TDChr[testDifficulty][TDChr[testDifficulty].length-1]-1] + chrIdOffset[chrIdOffset.length-1];
		System.out.println("Last chr and tot length=" +  Integer.toUnsignedString(lSeq));
		System.out.println("**End ChrIdOffset**" );
	}

	/** 
	 * getter for chrNA
	 * @return chrNA
	 */			
	public int[] getChrNA() {
		return chrNA;
	}

	/** 
	 * setter for chrNA
	 * @param chrNA - if NN blocks are loaded from a file
	 */					
	public void setChrNA(int[] cNA) {
		this.chrNA = cNA;
	}
	
	/** 
	 * setter for chrAA
	 * @param chrAA - if chromatid sequence in B is loaded from a file
	 */					
	public void setChrAA(byte[] cAA) {
		this.chrAA = cAA;
	}
	
	/** 
	 * getter for chrAA
	 * @return chrAA
	 */				
	public byte[] getBChrSequence() {
		return chrAA;
	}
	
	/** 
	 * getter for total Length 
	 * @return totLength
	 */					
	public int getTotalLength() {
		return totLength;
	}
	
	/** 
	 * Collect N using chrIdPtr and chrNA
	 *
	 * @param chrSeq - (input) An array with the chromatid sequence
	 * @param chrA - (output) A single long buffer is used for all the chrs. this is set in this method. (B rtpe)
	 * @param chrIdByteOffset - the byte offset into chrA where the the current chromatid is to be converted into B
	 *
	 * @return none
	 *         sets class variables chrNA and chrNCtr for NN Blocks
	 */
	void convChrToB(String[] chrSeq,byte[] chrA,int chrIdByteOffset) {
		
		int stB = chrIdByteOffset;   //byte array counter	
		int gCtr = chrIdByteOffset * 4; // # genome counter
		
		boolean startN = false;
		boolean foundN = false;
		
		BinaryRep br = new BinaryRep();
		
		for (int rc =0;rc < chrSeq.length; rc++) {
			int btCtr = 0;
			int limit = 0;

			for (int i =0; i< chrSeq[rc].length(); i+= 4) {
				int shift = 6;
				limit =  chrSeq[rc].length() - i;
				if (limit > 4)
					limit = 4;

				byte setB;
				for (int j=0; j< limit ; j++) {
					setB = 0b00000000;
					
					int c = (int)chrSeq[rc].charAt(i+j);
					
					if (c == 65) { //A
						setB = (byte) (g_vals[0]  << shift);
						foundN = false;
							
					} else if (c == 67) {
						setB = (byte) (g_vals[1]  << shift);
						foundN = false;
					}else if (c == 71) {
						setB = (byte) (g_vals[2]  << shift);
						foundN = false;
					}else if (c == 84) {
						setB = (byte) (g_vals[3]  << shift);
						foundN = false;
						
					}else {
						setB = (byte) (g_vals[0]  << shift);
						//set N
						foundN = true;
						if (!startN) {
							startN = true;
							chrNA[chrNCtr] = (gCtr + i + j);
							chrNCtr++;
						}
					}
					if (!foundN) {
						if (startN) {
							startN = false;
							chrNA[chrNCtr] = (gCtr + i + j );
							chrNA[0] += 2;
							chrNCtr++;
						}
					}
					chrA[stB+btCtr] |= setB;
					shift -= 2;
					
				}//endFor j
				btCtr += 1;
			}//endFor i
			stB += chrSeq[rc].length() / 4;

			gCtr += chrSeq[rc].length();
		}//endFor rc
		if (startN) {
			chrNA[chrNCtr] = gCtr;
			chrNA[0] += 2;
			chrNCtr++;
		}
		return;
	}

	/** 
	 * Initialization - sets up class variables for reading and creating the B file from the chromatid sequence
	 *                  and for alignment
	 *
	 * @param testDifficulty - 0 - small, 1 - medium, 2 - large
	 *                         testDifficulty 0 uses 1 chr, 1 uses 3 chromosomes, 2 uses ALL (24) chrs and 3 - uses 1 chr
	 *
	 * @return any valid integer
	 */
	public int initTest(int testDifficulty) {
		//System.out.println("DNAS:initTest() diff=" + testDifficulty);
		this.testDifficulty = testDifficulty;
		
		initReadBFile(testDifficulty);
		initCreateBFile(testDifficulty);
		initAlign(testDifficulty);
		return 0;
	}

	/** 
	 * Calculate total length required for chrAA - calc the length of all the chrs + GAP + max adjustment length of 4 bytes 
	 *
	 * @param testDifficulty - 0 - small, 1 - medium, 2 - large
	 *
	 * @return None
	 *         sets class variables chrIdA and totLength
	 */	
	void initReadBFile(int testDifficulty) {
		this.chrIdA = new int[TDChr[testDifficulty].length];

		this.totLength = (TDChr[testDifficulty].length -1) * (GAPBETCHRBYTE + 4);
		for (int x : TDChr[testDifficulty])
			this.totLength += chrLnSeq[x-1];
	}

	/** 
	 * Allocate chrAA and chrNA
	 *
	 * @param testDifficulty - 0 - small, 1 - medium, 2 - large
	 *
	 * @return None
	 *         sets class variables chrIdA and totLength
	 */		
	void initCreateBFile(int testDifficulty) {
		//could allocate individual length for each chr id
		//Currently max number of NN blocks is less than 170 (start,end) pairs per chromosome
		this.chrNA = new int[TDChr[testDifficulty].length * 340];

		this.chrAA = new byte[this.totLength];
	}

	/** 
	 * Initialization - if chromatid B sequence and chrNA are read from a preprocessed file
	 *
	 * @param testDifficulty - 0 - small, 1 - medium, 2 - large
	 *                         testDifficulty 0 uses 1 chr, 1 uses 3 chromosomes, 2 uses ALL (24) chrs and 3 - uses 1 chr
	 *
	 * @return none
	 */
	void initExt(int testDifficulty) {
		for (int i=0;i<TDChr[testDifficulty].length;i++) {
			chrIdA[chrIdPtr] = TDChr[testDifficulty][i];
			chrIdPtr++;
		}
	}
	
	/** 
	 * Initialization - if chromatid B sequence for all 24 chromosomes and chrNA are read from a preprocessed file
	 *                  and need to index a specific chromatid to test alignment against that chromatid
	 *
	 * @param testDifficulty - 0 - small, 1 - medium, 2 - large
	 *
	 * @return none
	 *         Class variable chrNAChr is updated here
	 */
	void initChrNAOffset(int testDifficulty) {
		chrNAChr = new int[(TDChr[testDifficulty].length *2) + 1];
		//set length to this val
		chrNAChr[0] = TDChr[testDifficulty].length *2;
		
		//traverse through chrNA and set the offsets
		int NStPos = 1;
		int currj = 1;
		for (int i=0;i<TDChr[testDifficulty].length;i++) {
			int chkVal = chrIdOffset[i];
			for (int j =currj;j<chrNA[0]; j+=2) {
				if (chrNA[j] == chkVal) {
					currj = j;
					chrNAChr[(i*2)+1] = currj;
					break;
				}
			}
			//update length;
			if (i > 0) {
				chrNAChr[((i-1)*2)+2] = currj - NStPos - 2; //subtract two for the filler in the GAP and onre for initial value
				NStPos = currj;
				currj += 2;
			}
		}
		//update last length
		chrNAChr[(TDChr[testDifficulty].length *2)] = chrNA[0] - NStPos;
	}
	
	/** 
	 * Initialization - variables needed for performing alignment
	 *                  Allocate probables 
	 *
	 * @param testDifficulty - 0 - small, 1 - medium, 2 - large
	 *
	 * @return none
	 *         Class variable chrIdOffset needed for alignment is set here
	 */
	void initAlign(int testDifficulty) {	
		chrIdOffset = new int[TDChr[testDifficulty].length];
		for (int i=1;i<TDChr[testDifficulty].length;i++) {
			int lSeq = chrLSeq[TDChr[testDifficulty][i-1]-1];
			chrIdOffset[i] = chrIdOffset[i-1] + lSeq + (4 - (lSeq % 4)) + GAPBETCHR;
		}
		
		for (int i=0; i< (splitA.length * 2); i++)
			for (int j=0; j< NUMPROBABLES; j++) {
				probablesM1[i][j] = new Probable();
				probablesM2[i][j] = new Probable();
			}
	}
	
	/** 
	 * Will be called 1, 3, or 24 times according to the test difficulty. Return can be any valid integer
	 *
	 * @param chromatidSequenceId - chromatid number
	 * @param chromatidSequence - nucleotide sequence of this chromatid
	 *
	 * @return any valid integer
	 *         Class variable chrIdA, chrIdPtr, chrAA are set here
	 */
	public int passReferenceGenome(int chromatidSequenceId, String[] chromatidSequence) {
		if (dbg > 1)
			System.out.println("DNAS:passReferencegenome chr=" + chromatidSequenceId + " chrIdPtr = " + chrIdPtr + " chrnCtr = " + chrNCtr );	
		//Add to chrIdA
		this.chrIdA[chrIdPtr] = chromatidSequenceId;
		
		//Convert to B and collect N Positions
		//If this is not the first chromosome, then insert NN Blocks
		if (chrIdPtr > 0) {
			chrNA[0] += 2;
			//start at the end of the previous chr
			chrNA[chrNCtr] = chrIdOffset[chrIdPtr-1] + chrLSeq[chrIdA[chrIdPtr - 1] - 1];;
			chrNCtr++;
			//end at the beginning of the new chr
			chrNA[chrNCtr] = chrIdOffset[chrIdPtr];
			chrNCtr++;
		}
		convChrToB(chromatidSequence,this.chrAA,Integer.divideUnsigned(chrIdOffset[chrIdPtr],4));
		
		this.chrIdPtr++;
		
		return 0;
	}

	/** 
	 * Any preprocessing before the reads are passed to performAlignment
	 * The chromatid is hashed here into a gigantic array
	 *
	 * @param none
	 * @return hashcount
	 */
	public int preProcessing() {
		return preProcessingExt(chrAA,chrNA);
	}

	/** 
	 * Call to Hashing if the chromatid is read from a preprocessed B file
	 * The chromatid is hashed here into a gigantic array
	 *
	 * @param chrAALocal - passed in chrAA (cromatid sequence in B read from a file
	 * @param chrNALocal - passed in chrNA (NN blockk array) read from a file
	 * @return hashcount
	 */
	public int preProcessingExt(byte[] chrAALocal,int[] chrNALocal) {
		return preProcessingExtChr(chrAALocal,chrNALocal,0);
	}

	/** 
	 * Call to Hashing if the chromatid is read from a preprocessed B file
	 * The chromatid is hashed here into a gigantic array
	 *
	 * @param chrAALocal - passed in chrAA (cromatid sequence in B read from a file
	 * @param chrNALocal - passed in chrNA (NN blockk array) read from a file
 	 * @param chrNo - a single chromatid id 
	 *                if 0 - all the chromatids loaded will to be hashed
	 * @return hashcount
	 */
	public int preProcessingExtChr(byte[] chrAALocal,int[] chrNALocal,int chrNo) {
		if (dbg > 1)
			System.out.println("DNASequencing:preProcesingExtChr() chrNo= " + chrNo);
		int hCtr = 0;
		hInd.init(this.k);
		int lastChr = 0; 
		int n = 0; //the end sequence number , if 0 the total length 
		int NNStart=0, NNLength = 0; //the offsets into the chrNAChr array for this chr ID

		//For all chr
		//total length incl insertions of N blocks between chromosomes
		if (chrNo == 0) {
			lastChr = chrIdA[chrIdA.length-1];
			n = chrIdOffset[chrIdOffset.length - 1]+chrLSeq[lastChr-1];
		}else {
			lastChr = chrIdA[chrNo-1];
			n = chrIdOffset[lastChr-1] + chrLSeq[lastChr-1];
			NNStart = chrNAChr[((lastChr-1)*2)+1];
			NNLength = chrNAChr[((lastChr-1)*2)+2];
			if (dbg > 1)
				System.out.println("DNASequencing:preProcesingExtChr() NNStart = " + NNStart + " NNLength =" + NNLength);
		}
		if (chrNo == 0)
			hCtr = hInd.rollingHash(chrAALocal,chrNALocal,n,this.k); 
		else {
			hCtr = hInd.rollingHashChr(chrAALocal,chrNALocal,chrIdOffset[lastChr-1],n,k,NNStart,NNStart + NNLength -1); 
		}
		return hCtr;
	}
	
	/**
	*   Will be called once per test size.
	*
    *   @param N is the number of reads in the test.
    *   Both array will contain exactly N elements as well as the return array. Elements at (2*i) and (2*i+1) are representing a read pair.
    *    
    *   @return result Array with alignment information
	*   Format of the return string:
    *
    *   ReadName, ChromatidSequenceId, Start position (1-based), End position (1-based), strand (reference [+] or reverse complement [-]), mapping confidence score
    *   sim1/1,1,10000,10149,+,0.89 
    *   sim1/2,1,10500,10652,-,0.82
    *   sim2/1,1,11500,11649,-,0.75
    *   sim2/2,1,11300,11449,+,0.63
	**/
	public String[] getAlignment(int N, double normA, double normS, String[] readName, String[] readSequence) {
		BinaryRep br = new BinaryRep();
		
		String[] res = new String[N];
		byte[] patA = new byte[READLENGTH];
		byte[] patB = new byte[READLENGTH];
		byte[] rpatA = new byte[READLENGTH];
		byte[] rpatB = new byte[READLENGTH];
		
		
		for(int i=0; i<N; i+=2) {
			//For debug
			/**
			if (i > 12) {
				res[i] =  readName[i] + ",20,26871291,26871440,-,0";
				res[i+1] =  readName[i+1] + ",20,26871291,26871440,-,0";
				continue;
			}
			**/
			//if (i != 176)
			//	continue;
				
			patA = br.convGLineToBBytes(readSequence[i],patA);
			patB = br.convGLineToBBytes(readSequence[i+1],patB);
			
			rpatA = br.reverseComplementBBytes(patA,rpatA);
			rpatB = br.reverseComplementBBytes(patB,rpatB);

			//if (dbg > 1)
				//System.out.println("i-" + i + " patA=" + BinaryRep.getString(patA,10) + " patB =" + BinaryRep.getString(patB,10)); 
			
			//Initialise probable counts
			for (int k=0; k< numProbablesM1.length; k++) {
				numProbablesM1[k] = 0;
				numProbablesM2[k] = 0;
			}
			
			//Get probables for patA,rpatB (patternA and reverse patternB)
			int numPairsM1 = matchPat(patA,rpatB,numProbablesM1,probablesM1,i,0);
			
			//Get probables for rpatA,patB
			int numPairsM2 = matchPat(rpatA,patB,numProbablesM2,probablesM2,i,1);
			//Deduce the probables for the pair based on patA and patB probables
			Probable prob = deduceProbable(numPairsM1,numPairsM2,i);
			
			//Convert to result Array with offsets calculated using chrIdOffset
			int tChr=0,treada=0,treadb = 0;
			boolean found = false;
			for (int x = 1;x<chrIdOffset.length;x++) 
				if (Integer.compareUnsigned(prob.reada,chrIdOffset[x]) < 0) {
					found = true;
					tChr = chrIdA[x-1];
					treada = prob.reada - chrIdOffset[x-1];
					treadb = prob.readb - chrIdOffset[x-1];
					break;
				}
			if (!found) {
				tChr = chrIdA[chrIdA.length - 1];
				treada = prob.reada - chrIdOffset[chrIdOffset.length - 1];
				treadb = prob.readb - chrIdOffset[chrIdOffset.length - 1];
			}
			prob.setRealChr(tChr,treada,treadb);
			
			String[] resStr = prob.getResStr();
			res[i] = readName[i] + resStr[0];
			res[i+1] = readName[i+1] + resStr[1];
			
			if (dbg > 0) {
				System.out.println(res[i]);
				System.out.println(res[i+1]);
			}

			//if (i > 0)
			//	break;
		}
		return res;
	}
	
	/**
	 *   call to free allocated memory
	 */
	public void exit() {
		hInd.free();
	}
	
	//calculate distance from ideal distance
	Probable getBest(int[] numProbables,Probable[][] probables,int ind) {
		float[] score = new float[numProbables[ind]];
		float sc;
		for (Probable p: probables[ind]) {
			//sc = (450 - p.diff) / (float)(450+1);
			sc = (IDEALREADGAP - p.diff) / (float)(IDEALREADGAP+1);
		}
		return probables[ind][0];
	}
	
	//Input: Probable resp.conf contains the distance score
	//                resp.countA and countB contains num matches
	Probable calcConf(Probable resp) {
		if (dbg > 0) {
			System.out.println("calcConf");
			resp.printDebug();
		}
		//Calc confidence here
		
		return resp;
	}

	//Selection
	Probable selectProbable(int[] numProbables,Probable[][] probables,int ind,int readNo) {
		Probable resp = probables[ind][0];
		if (numProbables[ind] == 1) {
			//get next probable
			//Check index of second place
			int nextInd = ind - 1; 
			while ((nextInd > 0) && (numProbablesM1[nextInd] == 0))
				nextInd--;
			
			if (dbg > 0)
				if (nextInd > 0)
					System.out.println("ind=" + ind + " nextInd=" + nextInd + " readNo=" + readNo);
					
			resp.conf = ((float)ind - nextInd) / (4 * (float)ind);
			resp.conf = resp.conf + (float)0.75;

		}if (numProbables[ind] > 1) {
			resp = getProbableDiff(numProbables,probables,ind,readNo);
			// If second place in same index, then conf  = 0.1
			resp.conf = (float)0.1;
		}
		return resp;
	}
	
	//get probables based on difference from IDEALREADGAP
	Probable getProbableDiff(int[] numProbables,Probable[][] probables,int ind,int readNo) {
		//All top go here
		int[] respA = new int[NUMPROBABLES];
		
		int numR = 0;
		
		float sc;// = 0.0;
		Probable resp = probables[ind][0];
		int diff = 10000; //some large number
		int currDiff = 0;
		Probable p = resp;
	 
		//Use distance forom 450 as the tie breaker
		for (int i=0; i< numProbables[ind];i++) {
			p = probables[ind][i];
			currDiff = Math.abs(IDEALREADGAP - p.diff);
			p.conf = (float)1.0 - (float)currDiff / READGAPDEN;
			if (dbg > 0) {
				p.printDebug();
				System.out.println("currDiff =" + currDiff);
			}
			
			if (currDiff <= diff) {
				resp = p;
				diff = currDiff;
				respA[numR] = i;
				numR++;
			}
		}
		
		//if more than one resp with equal diff, check popularity of strings
		if (numR > 1) {
			if (dbg > 1)
				System.out.println("numR= " + numR + " readNo=" + readNo);
			if (dbg > 1)
				for (int i=0; i< numR; i++)
					System.out.println("num= " + respA[i]);

			/** TO DO
			for (int i=0; i< numR; i++) {
				int readA = probables[ind][respA[i]].reada;
				int readB = probables[ind][respA[i]].readb;
				int totA = hInd.stringPopularity(chrAA[0],readA,150,splitA,45);
				int totB = hInd.stringPopularity(chrAA[0],readB,150,splitA,45);
				System.out.println("readA=" + readA + " totA=" + totA + " readB=" + readB + " totB=" + totB);
				resp = resp;
			}
			**/
		}
		return resp;
	}
	
	//selection based on fist match
	Probable getFirst(int[] numProbables,Probable[][] probables,int ind) {
		Probable p = probables[ind][0];
		if (numProbables[ind] > 1)
			p.conf = (float)0.1;
		else
			p.conf = (float)0.9;
		return p;
	}
	
	/**
	 * Given the list of probables, select the one based on First match
	 * uses class variables - numProbablesM1[] and numProbablesM2[] with the probables indexed by number of matches
	 *  
	 * @param numPairsM1 the number of pairs in M1 cycle ( not used)
	 * @param numPairsM2 the number of pairs in M2 cycle ( not used)
	 * @param readNo the read number - used only for printing DBG statement
	 * 
	 * @return the Probable match 
	 */	
	Probable deduceProbable(int numPairsM1,int numPairsM2,int readNo) {
		int type = 0;
		int indM1 = numProbablesM1.length - 1;
		while ((indM1 > 0) && (numProbablesM1[indM1] == 0))
			indM1--;

		int indM2 = numProbablesM2.length - 1;
		while ((indM2 > 0) && (numProbablesM2[indM2] == 0))
			indM2--;
		
		if (dbg > 1)
			System.out.println(":deduceProbable indM1 = " + indM1 + " indM2= " + indM2 + " readNo = " + readNo);
		if (dbg > 0)
			if (indM1 > 0 && indM2 > 0)
				System.out.println(":deduceProbable indM1 = " + indM1 + " indM2= " + indM2 + " readNo = " + readNo);
		
		Probable prob;
		if (indM1 > indM2) {
			type = 0;
			//prob = getProbableDiff(numProbablesM1,probablesM1,indM1);
			prob = selectProbable(numProbablesM1,probablesM1,indM1,readNo);
			
		}else if (indM1 < indM2) {
			type = 1;
			//prob = getProbableDiff(numProbablesM2,probablesM2,indM2);
			prob = selectProbable(numProbablesM2,probablesM2,indM2,readNo);
		} else {
			type = 0;
			if (dbg > -1)
				System.out.println("DBG deduceProbable() indM1 == indM2 " + indM1 + " readNo = " + readNo);
			//Selecting M1 arbitrarily
			Probable p = probablesM1[indM1][0];
			p.conf = (float)0.0;
			return p;
			
		}
		prob.type = type;
		//TO DO calc conf
		//prob = calcConf(prob);		
		return prob;
	}

	/**
	 * Given the list of probables, select the one based on First match
	 * uses class variables - numProbablesM1[] and numProbablesM2[] with the probables indexed by number of matches
	 *  
	 * @param ppA pattern A
	 * @param ppB pattern B
	 * @param readNo the read number - used only for printing DBG statement
	 * 
	 * @return numPairs number of pairsthat meet criteria
	 *         Class variables probables, numProbables are updated here
	 */		
	int matchPat(byte[] ppA,byte[] ppB,int[] numProbables,Probable[][] probables,int readNo,int type) {
		//Get Probables for each subset
		if (dbg > 1) {
			System.out.println("probables ppA=");
			BinaryRep.printString(ppA);
			System.out.println("probables ppA readNo =" + readNo );
		}
		Map<Integer,int[]> probA = hInd.querySplit(ppA,splitA,this.k);
		compressProbables(probA,splitA.length);
		//HashIndex.printCounter(probA);

		if (dbg > 1) {
			System.out.println("probables ppB=");
			BinaryRep.printString(ppB);
			System.out.println("probables ppB readNo =" + readNo );
		}
		Map<Integer,int[]> probB = hInd.querySplit(ppB,splitA,this.k);
		compressProbables(probB,splitA.length);
		//HashIndex.printCounter(probB);
		
		//Get Freq Distb
		int[] freqProbA = new int[splitA.length];
		
		for (int[] x: probA.values()) 
			freqProbA[x[0]-1]++;
		int[] freqProbB = new int[splitA.length];
		for (int[] x: probB.values()) 
			freqProbB[x[0]-1]++;

		//HashIndex.printFreqProb(freqProbA);
		//HashIndex.printFreqProb(freqProbB);
		
		int[][] chrProbA = new int[splitA.length][];
		int[][] chrProbB = new int[splitA.length][];
		for (int i=0;i<splitA.length;i++) {
			chrProbA[i] = new int[freqProbA[i]+1];
			chrProbB[i] = new int[freqProbB[i]+1];
		}
		//Dbg print
		/**
		System.out.println("charProbA");
		for (int i=0; i< chrProbA.length;i++) 
			System.out.print(chrProbA[i].length + " ");
		System.out.println("chrProbB");
		for (int i=0; i< chrProbB.length;i++) 
			System.out.print(chrProbB[i].length + " ");
		System.out.println();
		**/
		
		//Collect values
		for (int key : probA.keySet()) {
			int val = probA.get(key)[0] - 1;
			chrProbA[val][0]++;
			chrProbA[val][chrProbA[val][0]] = key;
		}
		//Collect values
		for (int key : probB.keySet()) {
			int val = probB.get(key)[0] - 1;
			chrProbB[val][0]++;
			chrProbB[val][chrProbB[val][0]] = key;
		}
		
		//Dbg print
		/**
		if (dbg > 0) {
			System.out.println("charProbA");
			for (int i=0; i< chrProbA.length;i++) 
				System.out.print(chrProbA[i][0] + " ");
			System.out.println("chrProbB");
			for (int i=0; i< chrProbB.length;i++) 
				System.out.print(chrProbB[i][0] + " ");
			System.out.println();
		}
		**/
			
		//Collect pairs
		int iA=splitA.length-1;
		while (iA > 0 && freqProbA[iA] == 0)
			iA--;
		int iB=splitA.length-1;
		while (iB > 0 && freqProbB[iB] == 0)
			iB--;
		int iAInit = iA;
		int iBInit = iB;

		int numPairs=0;
		while (iA>=0) {
			if (dbg > 1)
				System.out.println("iA=" + iA + " iB="+iB);
			while (iB>=0) {
				if (dbg > 1)
					System.out.println("iB=" +iB);
	
				int a=1,b=1;
		
				while (a < chrProbA[iA].length && b < chrProbB[iB].length) {
					int itemA = chrProbA[iA][a];
					int itemB = chrProbB[iB][b];
					int diff =  Math.abs(itemA - itemB);

					if (diff > MINREADGAP && diff < MAXREADGAP) {
						if (dbg > 0)
							System.out.println("diff=" + diff + " reada=" + itemA + " readB=" + itemB + " iA=" + iA + " iB=" + iB + " tot=" + (iA+iB+2));
						int tot = iA+iB;

						//Collect only 10 probables in any freq
						if (numProbables[tot] < NUMPROBABLES) {

							probables[tot][numProbables[tot]].setValues(diff,itemA,itemB,iA,iB,tot,type);
							numProbables[tot] += 1;
							
							//Increment number of probables found
							numPairs += 1;
							//Dbg
							if (dbg > 0)
								if (numProbables[tot] == NUMPROBABLES)
									System.out.println("DBG numprobables[tot] is max= " + tot + " readNo = " + readNo);
						}else
							break;
					}
					if (Integer.compareUnsigned(itemA,itemB) < 0)
						a++;
					else
						b++;
		
				}//endfor a &&b
			iB--;
			}//endWhile iB
			iA--;
			iB= iBInit;
						
		}//endWhile iA and iB
		return numPairs;
	}
	
	/**
	 * if the reads in the probables are less than tolerance apart, compress them into a single read and increase the count for that read
	 *  
	 * @param prob (read,count) pair
	 * @param maxVal count cannot exceed maxVal (split.length is the maxVal)
	 * 
	 * @return none
	 *         Map prob is updated here (reads are merged and count is increased)
	 */		
	void compressProbables(Map<Integer,int[]> prob,int maxVal) {
		int tolerance = 3;
		int prevKey = -1000;
		ArrayList<Integer> remList = new ArrayList<Integer> ();
		
		for (int x : prob.keySet()) {
			if (Integer.compareUnsigned((x - prevKey),tolerance) < 0) {
				//add ctr of current key to prevKey
				int[] val= prob.get(prevKey);
				int currVal = prob.get(x)[0];
				val[0] += prob.get(x)[0];
				if (val[0] > maxVal) 
					val[0] = maxVal;
				
				//Remove curr Key
				remList.add(x);
			}
			prevKey = x;
		}
		//Remove keys
		for (int r: remList)
			prob.remove(r);
	}

	
}