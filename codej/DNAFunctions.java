import java.time.Instant;
import java.time.Duration;
import java.time.temporal.Temporal;
import java.time.temporal.ChronoUnit;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

//For minisam file
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;

import java.io.BufferedWriter;
import java.io.FileWriter;

import java.util.Collections;


public class DNAFunctions {

	/**
	* Constants from the problem statement
	*/
	static final int MAX_POSITION_DIST = 300;
	static final String resDirPath = "../data/results/";
	static final String resExt = ".minis";


 	public void createBFile(int testDifficulty) {
		System.out.println("DNAFunctions:createBFile(): Start option=" + testDifficulty);
		if ((testDifficulty > 2) || (testDifficulty <0))
			return;
		//output where the B file will be stored
		String bPath = "../data/Bchroma/";
		String chrPath = "../data/chromatids/chromatid";
		String ext = ".b";
		
		int[] chr_ids = new int[1];	
		if(testDifficulty==0) {
			chr_ids = new int[]{20};
		} else if(testDifficulty==1) {
			chr_ids = new int[]{1,11,20};
		} else if(testDifficulty==2) {
			chr_ids = new int[24];
			for(int i=1; i<=24; ++i) 
				chr_ids[i-1] = i;		
		}
	
		Instant t11,t22;

		DNASequencing dnas = new DNASequencing();
		dnas.initTest(testDifficulty);
		dnas.printChrIdOffset(testDifficulty);

		// load chromatid	
		DataIO d = new DataIO();
		t11 = Instant.now();
		for(int chromatid_seq_id: chr_ids) {
			String[] chromatid_seq;
			String fname =  chrPath + String.valueOf(chromatid_seq_id) + ".fa";
			
			chromatid_seq = d.readFastaFile(fname);

			dnas.passReferenceGenome(chromatid_seq_id, chromatid_seq);	
			//System.out.println("createBFile(): chrID = " + chromatid_seq_id);
			//dnas.printChrNA(0);
		}
		t22 = Instant.now();
		System.out.println("total time for passReferenceGenome ="  + Duration.between(t11, t22).toString());		
		
		int[] chrNA = dnas.getChrNA();
		//dnas.printChrNA(0);
		String chrNAPath = bPath + "Chr-NN" + testDifficulty + ext;
		System.out.println("DNAFunctions:createBFile(): chrNAPath = " + chrNAPath);
		//d.saveIntFile(chrNAPath,chrNA);
		
		byte[] chrSeq = dnas.getBChrSequence();
		String chrSeqPath = bPath + "Chr-Seq" + testDifficulty + ext;
		System.out.println("DNAFunctions:createBFile(): chrSeqPath = " + chrSeqPath);
		//d.saveBinaryFile(chrSeqPath,chrSeq);
		
	}
	
 	public void readBFile(int testDifficulty) {
		System.out.println("DNAFunctions:readBFile() Start option=" + testDifficulty);
		DNASequencing dnas = new DNASequencing();
		
		readBFileIndex(testDifficulty,dnas);
		//dnas.printChrNA(0);
	}

 	public void readBFileIndex(int testDifficulty,DNASequencing dnas) {
		System.out.println("DNAFunctions:readBFileIndex(): Start option=" + testDifficulty);
		readBFileIndexChr(testDifficulty,dnas,0);
	}
	
	
	//if chrNo = 0, then all chrs
 	public void readBFileIndexChr(int testDifficulty,DNASequencing dnas,int chrNo) {
		System.out.println("DNAFunctions:readBFileIndexChr(): Start option=" + testDifficulty + " chrNo= " + chrNo);
		//output where the B file will be stored
		String bPath = "../data/Bchroma/";
		String chrPath = "../data/chromatids/chromatid";
		String ext = ".b";
		
		int[] chr_ids = new int[1];	
		if(testDifficulty==0) {
			chr_ids = new int[]{20};
		} else if(testDifficulty==1) {
			chr_ids = new int[]{1,11,20};
		} else if(testDifficulty==2) {
			chr_ids = new int[24];
			for(int i=1; i<=24; ++i) 
				chr_ids[i-1] = i;		
		}
		Instant t1,t2;
		dnas.initReadBFile(testDifficulty);
		
		//Read NN file
		DataIO d = new DataIO();
		String chrNAPath = bPath + "Chr-NN" + testDifficulty + ext;
		int[] NNln = new int[1];
		NNln = d.loadFromIntFile(chrNAPath,0,1);
		int ln = NNln[0] + 1;
		int[] chrNA = d.loadFromIntFile(chrNAPath,0,ln);	
		dnas.setChrNA(chrNA);
		//dnas.printChrNA(0,0);

		int totLength = dnas.getTotalLength();
		String chrSeqPath = bPath + "Chr-Seq" + testDifficulty + ext;		
		byte[] chrAA = d.loadBinaryFile(chrSeqPath,0,totLength);

		dnas.initAlign(testDifficulty);
		dnas.initExt(testDifficulty);
		
		dnas.initChrNAOffset(testDifficulty);
		//dnas.printChrIdOffset(testDifficulty);
		//dnas.printChrNAChr();

		System.out.println("DNAFunctions:readBFile() Start preProcessing");
		
		int hCtr = 0;
		t1 = Instant.now();
		hCtr=dnas.preProcessingExtChr(chrAA,chrNA,chrNo);
		t2 = Instant.now();
		System.out.println("readBFile() time for preProcessing hashcount=" + hCtr + "time=" + Duration.between(t1, t2).toString());		
		
		
	}
	
	void testCreateNNChr(int testDifficulty) {
		System.out.println("DNAFunctions:testCreateNNChr(): Start");
		//output where the B file will be stored
		String bPath = "../data/Bchroma/";
		String chrPath = "../data/chromatids/chromatid";
		String ext = ".b";

		DNASequencing dnas = new DNASequencing();
	
		dnas.initReadBFile(testDifficulty);
		
		//Read NN file
		DataIO d = new DataIO();
		String chrNAPath = bPath + "Chr-NN" + testDifficulty + ext;
		int[] NNln = new int[1];
		NNln = d.loadFromIntFile(chrNAPath,0,1);
		int ln = NNln[0] + 1;
		int[] chrNA = d.loadFromIntFile(chrNAPath,0,ln);	
		dnas.setChrNA(chrNA);
		//dnas.printChrNA(NNStart,NNLength);

		dnas.initAlign(testDifficulty);
		dnas.initExt(testDifficulty);
		
		dnas.initChrNAOffset(testDifficulty);
		//dnas.printChrIdOffset(testDifficulty);
		dnas.printChrNAChr();
		
		
	}

	public void perfAlignOne(int testDifficulty, int readOption) {
		System.out.println("DNAFunctions:perfAlignOne(): Start Chr option=" + testDifficulty + " Readfile option = " + readOption);
		DNASequencing dnas = new DNASequencing();
		//read che file and index it
		readBFileIndex(testDifficulty,dnas);
		//dnas.printChrNA(0);
		
		perfAlign(testDifficulty,readOption,dnas);
		dnas.exit();
	}
	
	public void perfAlign(int testDifficulty, int readOption,DNASequencing dnas) {
		System.out.println("DNAFunctions:perfAlign(): Start Chr option=" + testDifficulty + " Readfile option = " + readOption);
		
		String[] readFiles = {"../data/exampletest/small5.fa1",
		"../data/exampletest/small5.fa2",
		"../data/exampletest/medium5.fa1",
		"../data/exampletest/medium5.fa2"};
		
		String fa1Path= readFiles[readOption * 2];
		String fa2Path= readFiles[(readOption * 2) + 1];
		
		Instant t1,t2;

		
		//read align file
		// load reads
		String resFile = formResFilename(testDifficulty,readOption);		
		
		String[] readName = new String[2];
		String[] readSequence = new String[2];
		try {
			//reading file line by line
			FileInputStream fin1 = new FileInputStream(fa1Path);
			BufferedReader reader1 = new BufferedReader(new InputStreamReader(fin1));

			FileInputStream fin2 = new FileInputStream(fa2Path);
			BufferedReader reader2 = new BufferedReader(new InputStreamReader(fin2));
			
			//write result file
			FileWriter fout = new FileWriter(resFile);
			BufferedWriter bwtr = new BufferedWriter(fout);
			String content = "This is the content to write into file\n";
			
			String line1 = reader1.readLine(); //name line
			String line2 = reader2.readLine(); //name line
			
			t1 = Instant.now();
			while(line1 != null && line2 != null){
				//Add name
				readName[0] = line1.substring(1,line1.length());
				readName[1] = line2.substring(1,line2.length());

				line1 = reader1.readLine(); //read seq
				line2 = reader2.readLine();
				readSequence[0] = line1;
				readSequence[1] = line2;
				
				//Call align			
				int nreads = readName.length;
				// compute alignments
				
				String[] results = dnas.getAlignment(nreads, 0.5, 0.5, readName, readSequence);
				//Save result
				for (int i=0;i<2;i++) {
					//System.out.println("res=" + results[i]);
					bwtr.write(results[i]);
					bwtr.newLine();
				}
				
				line1 = reader1.readLine(); //read name 
				line2 = reader2.readLine();
			}           
			t2 = Instant.now();
			System.out.println("time for getAlignment readname chrOpt =" + testDifficulty + " read option = " + readOption + " time=" + Duration.between(t1, t2).toString());		

			reader1.close();
			reader2.close();
			fin1.close();
			fin2.close();
			bwtr.close();
			fout.close();
			
			
		}catch(Exception e) {
			System.out.println(e);
			e.printStackTrace();
		}
		
		
	}
	
	String formResFilename(int testDifficulty,int readOption) {
		String resFile = resDirPath + "res-c" + testDifficulty + "-" + readOption + resExt;
		return resFile;
	}

	
	public void compareResults(int testDifficulty,int readOption) {
		String[] truthFiles = {"../data/exampletest/small5.minisam",
		"../data/exampletest/medium5.minisam",
		"../data/exampletest/large5.minisam"};
		
		String minisamPath = truthFiles[readOption];
		
		String resFile = formResFilename(testDifficulty,readOption);


		Map<String, Position> truth = parse_truth(minisamPath);	
		System.out.println("Truth size=" + truth.size());
		//read res file
		DataIO d = new DataIO();
		String[] results = d.readFile(resFile);
		List<ReadResult> readResults = build_read_results(truth, results);
		System.out.println("read results size=" + readResults.size());
	}
	
	
	/**
	* Read a minisam file and build a map of ground truth
	* @param path	the path of the minisam file storing the ground truth 
	* @return a map[read_name] = read_Position
	*/
	Map<String, Position> parse_truth(final String path) {
		Map<String, Position> res = new HashMap<String,Position>();

		try {
			FileInputStream fin = new FileInputStream(path);
			BufferedReader reader = new BufferedReader(new InputStreamReader(fin));
			
			String line = reader.readLine(); 
			int tno = 0;
			char prevSt = ' ';
			while(line != null){
				String[] tokens = line.split(",");
				Position posn = new Position(Integer.parseInt(tokens[1]),Integer.parseInt(tokens[2]),Integer.parseInt(tokens[3]),tokens[4].charAt(0));
				if ((tno % 2) == 0 )
					prevSt = tokens[4].charAt(0);
				else 
					if (prevSt == tokens[4].charAt(0))
						System.out.println("Two strands same dirn at posn=" + tno);
				
				tno += 1;
				
				res.put(tokens[0],posn);
				line = reader.readLine();
			}           
			reader.close();
			fin.close();
			return res;

		}catch(Exception e) {
			System.out.println(e);
			return null;
		}
	}

/**
 * For each string of the results vector, build a read result {confidence, r}
 * @param truth		the map of ground truth position for each read
 * @param results	the vector of results as return by getAlignment
 * @return a vector of ReadResult, that is {confidence, r}
 */
	List<ReadResult> build_read_results(final Map<String, Position> truth, final String[] results) {
		List<ReadResult> read_results = new ArrayList<ReadResult> ();
		int n = results.length;
		int correct = 0;
		
		for(int i=0; i<n; ++i) {
			String[] tokens = results[i].split(",");
			Position posn = truth.get(tokens[0]);
			//System.out.println("posn key=" + tokens[0] + " chrName=" + posn.chrName + " from=" + posn.from + " to=" + posn.to + " strand=" + posn.strand);

			int r = 1;
			r = (Integer.parseInt(tokens[1]) == posn.chrName)? r : 0;
			r = (tokens[4].charAt(0) == posn.strand)? r : 0;
			int start0 = Integer.parseInt(tokens[2]);
			int start1 = posn.from;
			
			r = (Math.abs(start0-start1) < MAX_POSITION_DIST ) ? r : 0;
			double confidence = Double.parseDouble(tokens[5]);
			read_results.add(new ReadResult(confidence,r));
			correct += r;
		}
		System.out.println("Number of correct answers: " + correct + "/" + n + " = " + (double)correct/(double)n);
		return read_results;
	}	
	
	void customTestOne() {
		System.out.println("DNAFunctions:customTest(): Start Using chr 1 and reads 0 and 1");
		DNASequencing dnas = new DNASequencing();
		//read chr file and index it
		readBFileIndex(1,dnas);
		
		//test with chrOption 1 and files - 0 and 1
		perfAlign(1,0,dnas);
		perfAlign(1,1,dnas);

		dnas.exit();

		
		compareResults(1,0);
		compareResults(1,1);
	}
	
	void customTest24() {
		System.out.println("DNAFunctions:customTest24(): Start Using chr ALL and reads 0 and 1");
		DNASequencing dnas = new DNASequencing();
		//read chr file and index it
		readBFileIndex(2,dnas);
		
		//test with chrOption 2 and files - 0 and 1
		perfAlign(2,0,dnas);
		perfAlign(2,1,dnas);

		dnas.exit();
		
		compareResults(2,0);
		compareResults(2,1);
	}

	void customTestChr20() {
		System.out.println("DNAFunctions:customTestChr20(): Start Using chr 20 and reads 0");
		DNASequencing dnas = new DNASequencing();
		//read chr file and index it
		readBFileIndexChr(2,dnas,20);
		
		//test with chrOption 2 and files - 0 
		perfAlign(2,0,dnas);

		dnas.exit();
		
		compareResults(2,0);
	}
	

	/**
	* Main function: read the data, perform the DNA alignments and score results
	*                - An illustration of all the functions available here
	*                  1. Read HG chromosome file(s) and create B files
	*                     A B file is one that is encoded four DNA characters to a byte
	*/
    public static void main(String args[]) {
		/**
		 * options are 0,1,2 representing files (20, (1,11,20) ans (all 24) chr files
		 */
		int chrOption = 2;
		/**
		 * options are 0,1,2 representing with 19116,1878654 and LARGE number of paired end reads
		 */		
		int readOption =0;
		// perform test
		DNAFunctions d = new DNAFunctions();
		//d.createBFile(chrOption);
		
		//d.readBFile(chrOption);
		
		//d.perfAlignOne(chrOption,readOption);
		//d.compareResults(chrOption,readOption);
		//d.customTestOne();
		//d.customTest24();
		//d.customTestChr20();
		d.testCreateNNChr(2);

	return;
	}
}