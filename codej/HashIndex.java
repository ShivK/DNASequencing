//import java.util.HashMap;
import java.util.TreeMap;
import java.util.Map;
import java.util.stream.IntStream;

import java.lang.reflect.Field;

//Uses unsafe

/**
#B - the base of the numeral system, (B >= n)
# n the size of the alphabet
# and M - a big enough prime number
#Prime  179,424,673, 982,451,653, 373,587,883 (20,000,000),15,485,863
# 4 bytes max int = 4294967295
# 3 bytes maxint - 16777215
# 3.5 bytes - 268435455
# 3.25 bytes - 67108863
**/

public class HashIndex {
	//Hashing constants
    //int B = 4;
    //int M = 13;
	private static final int dbg = 0;

	private static final int B = 5;
	//int M = 15485863;
	private static final int M = 300000497; //300000007
	//private static final int M = 200000039;
	
	//private static final int M = 400000009;
	
	//This value should be less than 1000, the separation between chromosomes
	private static final int NNLimit = 4; //greater than NNLimit is considered as a N Block
	//private static final int NNLimit = 900; //greater than NNLimit is considered as a N Block
	
	static int hashCount; //unsigned int
	static int firstMapPosn;
	private static int[] mapStartPosn;
	
	private static sun.misc.Unsafe unsafe = getUnsafe();
	private static long mapValues;
	
	private static BinaryRep br;

    private static sun.misc.Unsafe getUnsafe() {
		try {
			Field f = sun.misc.Unsafe.class.getDeclaredField("theUnsafe");
			f.setAccessible(true);
			return (sun.misc.Unsafe) f.get(null);
		} catch (Exception e) {
			System.err.println("ERROR: getting unsafe: " + e.getMessage());
			return null;
		}
	}
	
	public static void init(int k) {
		//Not using k here
		br = new BinaryRep();
	}

	public static void free() {
		unsafe.freeMemory(mapValues);
	}
	
	private static int int_mod(int a,int b) {
		return (a % b + b) % b;
	}


	//m = len of kmer
    public static int expHash(byte[] pat,int startPosn,int m) {
        //# calculate the hash value of the pattern
        int hp = 0;
        for (int i=startPosn;i<startPosn+m; i++)
            hp = int_mod(hp * B + pat[i], M);
        return hp;
	}
	
	static void freqDistb(short[] mapT) {
		int[] freqA = new int[256];
		int num255 = 256;
		int num1024 = 1024;
		int num16 = 65536;
		int maxN = 0;;
		int f0=0,f255 = 0,f1024=0,f16=0;
		
		for (int i: mapT) {
			if (i < num255) {
				f0++;
				freqA[i]++;
			}
			if (i >= num255)
				f255++;
			if (i >= num1024)
				f1024++;
			if (i >= num16)
				f16++;
			if (i > maxN)
				maxN = i;
		}
				
		System.out.println("f0 = " + f0 + " f255= " + f255 + " f1024 = " + f1024 +  " f16 = " + f16); 
		System.out.println("maxN = " + maxN);
		
		for (int i=0;i<freqA.length;i++)
			if (freqA[i] !=0)
				System.out.println("i =" + i + " val =" +  freqA[i]);
		
	}
	

	public static short[] rollingHashInt1(byte[]tA,int[] NPos,int seqStartPosition,int n,int m,short[] mapT,int NNStart,int NNEnd) {

		//Start N Loop
		//Hash all subsequences that are separated by NN blocks
		int seqStPos = seqStartPosition; //unsigned
		int seqEndPos = 0; //unsigned
		int NIndex = NNStart - 1;
		//get Valid start of NN Block
		while (NIndex < NNEnd) { //NIndex < number of NN blocks in this chr
			//get valid NN Block
			if (Integer.compareUnsigned((NPos[NIndex+2] - NPos[NIndex+1]), NNLimit) > 0) {
				int stNN = NPos[NIndex+1];
				if (Integer.compareUnsigned(seqStPos, stNN) < 0) {
					seqEndPos = NPos[NIndex+1];
					if (!(Integer.compareUnsigned((seqEndPos - seqStPos),m) < 0)) {
						mapT =  rollingHashB1(tA,seqStPos,seqEndPos,m,mapT);
					}
				}
				seqStPos = NPos[NIndex+2];
			}	
			NIndex += 2;
		}

		//past all NN Blocks
		seqEndPos = n; //get the total legth of the chtomosome
		if (!(Integer.compareUnsigned((seqEndPos - seqStPos),m) < 0)) {
			mapT =  rollingHashB1(tA,seqStPos,seqEndPos,m,mapT);
		}
		//System.out.println("After first pass hashCount="+Integer.toUnsignedString(hashCount));
		
		//Get Freq Distb
		//freqDistb(mapT);
		
		mapStartPosn = new int[M];
		initMapStartPosn(mapT);
		//printMap(mapStartPosn);

		return mapT;
	}

	/**
      * int m : kmer
      * (unsigned) int n : Chr Length (not the number of bytes of the compressed tA) 
	  */
	public static int rollingHash(byte[]tA,int[] NPos,int n,int m) {
		return rollingHashChr(tA,NPos,0,n,m,1,NPos[0]);
	}

	//seqStartPosition
	// n - seqEndPosition
	public static int rollingHashChr(byte[]tA,int[] NPos,int seqStartPosition,int n,int m,int NNStart,int NNEnd) {
		if (dbg > 1) {
			System.out.println("HashIndex:rollingHash() N length=" + NPos[0] + "chr St posn = " + Integer.toUnsignedString(seqStartPosition) + " chr (end posn) length=" + Integer.toUnsignedString(n) + " kmer=" + m + " NNdiff=" + NNLimit + " M= " + M);
			System.out.println("HashIndex:rollingHash() NNStart=" + NNStart + " NN End=" + NNEnd);
		}
		short[] mapT = new short[M];
		mapT = rollingHashInt1(tA,NPos,seqStartPosition,n,m,mapT,NNStart,NNEnd);

		//Start N Loop
		//Hash all subsequences that are separated by NN blocks
		int seqStPos = seqStartPosition; //unsigned
		int seqEndPos = 0; //unsigned
		int NIndex = NNStart -1;		
		//Got all hashCounts, alloc memory (using unsafe)
		long lhc = Integer.toUnsignedLong(hashCount);
		lhc = lhc * Integer.BYTES;
		mapValues = getUnsafe().allocateMemory(lhc);
		
		//Fill up map positions
		//Start N Loop
		//Hash all subsequences that are separated by NN blocks
		//get Valid start of NN Block
		while (NIndex < NNEnd) {
			//get valid NN Block
			if (Integer.compareUnsigned((NPos[NIndex+2] - NPos[NIndex+1]), NNLimit) > 0) {
				int stNN = NPos[NIndex+1];
				if (Integer.compareUnsigned(seqStPos, stNN) < 0) {
					seqEndPos = NPos[NIndex+1];
					if (!(Integer.compareUnsigned((seqEndPos - seqStPos),m) < 0)) {
						mapT =  rollingHashB2(tA,seqStPos,seqEndPos,m,mapT);
					}
				}
				seqStPos = NPos[NIndex+2];
			}
			NIndex += 2;
		}
		//past all NN Blocks
		seqEndPos = n; //get the total legth of the chromosome
		if (!(Integer.compareUnsigned((seqEndPos - seqStPos),m) < 0)) {
			mapT =  rollingHashB2(tA,seqStPos,n,m,mapT);
		}
		
		
		/**
		System.out.println("MapValues=");
		printMap(mapValues);
		System.out.println("MapTotals=");
		printMap(mapT);
		System.out.println("MapStartPosn=");
		printMap(mapStartPosn);
		**/
		
		return hashCount;
	}
	
	/**
	 * Updates Class variable - firstMapPosn and mapStartPosn
	 */
	private static void initMapStartPosn(short[] mapT) {
		int startPosn = 0;
		//Update firstStartPosn
		for (int i=0;i<mapT.length; i++) 
			if (mapT[i] > 0) {
				mapStartPosn[i] = startPosn;
				startPosn += mapT[i];
				mapT[i] = 0;
				firstMapPosn = i;
				break;
			}
		//init rest
		for (int i=firstMapPosn+1;i<mapT.length; i++) 
			if (mapT[i] > 0) {
				mapStartPosn[i] = startPosn;
				startPosn += mapT[i];
				mapT[i] = 0;
			}
	}
	
	//Assume seqEndPos - seqStPos < m 
	/**
	  * (unsigned) int seqStPos
	  * (unsigned) int seqEndPos
	  */
    private static short[] rollingHashB1(byte[]chrA,int seqStPos,int seqEndPos,int m,short[] map) {
		
		byte[] tA = br.getBytesFromBB(chrA,seqStPos,m);
		int circBufPtr = 0;
		
        //calculate the hash value of the first segment 
        //of the text of length m
        int ht = 0;
		for (int i=0;i<m; i++)
            ht = int_mod(ht * B + tA[i], M);
		
		hashCount += 1;
		
		//Count
		map[ht] += 1;
		
		//if the diff id equal to m, no rolling hash needed
		if (Integer.compareUnsigned((seqEndPos - seqStPos),m) == 0)
			return map;

        //#start the "rolling hash" - for every next character in
        //#  the text calculate the hash value of the new segment
        //#  of length m; E = (Bm-1) modulo M

		//get the start byte point
		int bInd = Integer.divideUnsigned((seqStPos+m),4);
		//get the ith char in the mask
		int nChar = Integer.remainderUnsigned((seqStPos+m),4);
		byte rb;
		//get the byte
		byte b = chrA[bInd];

		//Calc E
        int ee = B % M;
        int E = ee;
		int ti = seqStPos;
		while (Integer.compareUnsigned(ti,(seqStPos+m-2)) < 0) {
			E = Integer.remainderUnsigned((E * ee),M);
			ti++;
		}

		for (int i=seqStPos+m; (Integer.compareUnsigned(i,seqEndPos) < 0);i++) {
			ht = int_mod(ht - int_mod(tA[circBufPtr] * E, M), M); //first char (i-m) is where the circBufPtr points
			
            ht = int_mod(ht * B, M);
			//Get next char
			rb = (byte)(((chrA[bInd] & br.g_masks[nChar])&0xff) >> ((3 - nChar) * 2));
			tA[circBufPtr] = rb;
			
			ht = int_mod(ht + tA[circBufPtr], M); //tA[circBufPtr] now contains the latest char
	
			//Inc byte counters
			circBufPtr = (circBufPtr + 1) % m;
			nChar += 1;
			if (nChar >= 4) {
				bInd += 1;
				nChar = nChar % 4;
				b = chrA[bInd];
			}
			map[ht] +=1;
			hashCount += 1;
		}
		return map;
	}
	
	//HashCount=63941582 hCtr =63941582
	/**
	  * Fill up hashpositions
	  * (unsigned) int seqStPos
	  * (unsigned) int seqEndPos
	  */
	private static short[] rollingHashB2(byte[]chrA,int seqStPos,int seqEndPos,int m,short[] map) {
		
		byte[] tA = br.getBytesFromBB(chrA,seqStPos,m);
		int circBufPtr = 0;

        //calculate the hash value of the first segment 
        //of the text of length m
        int ht = 0;
		for (int i=0;i<m; i++) 
            ht = int_mod(ht * B + tA[i], M);
		//hashCount += 1;
		
		//Count
		long offset = Integer.toUnsignedLong(mapStartPosn[ht] + map[ht]);
		unsafe.putInt(mapValues + offset * Integer.BYTES, seqStPos);
		
		map[ht] += 1;
		
		//if the diff id equal to m, no rolling hash needed
		if (Integer.compareUnsigned((seqEndPos - seqStPos),m) == 0)
			return map;

        //#start the "rolling hash" - for every next character in
        //#  the text calculate the hash value of the new segment
        //#  of length m; E = (Bm-1) modulo M
		//get the start byte point
		int bInd = Integer.divideUnsigned((seqStPos+m),4);
		//get the ith char in the mask
		int nChar = Integer.remainderUnsigned((seqStPos+m),4);
		
		//Calc E
        int ee = B % M;
        int E = ee;
		for (int i=seqStPos;(Integer.compareUnsigned(i,(seqStPos+m-2))<0);i++)			
            E = (E * ee) % M;
         
		for (int i=seqStPos+m; (Integer.compareUnsigned(i,seqEndPos) < 0);i++) {

			ht = int_mod(ht - int_mod(tA[circBufPtr] * E, M), M); //first char (i-m) is where the circBufPtr points
			
            ht = int_mod(ht * B, M);
			//Get next char
			tA[circBufPtr] = (byte)(((chrA[bInd] & br.g_masks[nChar])&0xff) >> ((3 - nChar) * 2));
			
			ht = int_mod(ht + tA[circBufPtr], M); //tA[circBufPtr] now contains the latest char

			//mapValues[mapStartPosn[ht] + map[ht]] = i-m+1; 
			long offs = Integer.toUnsignedLong(mapStartPosn[ht] + map[ht]);		
			
			unsafe.putInt(mapValues + offs * Integer.BYTES, i-m+1);

			map[ht] +=1;
			//hashCount += 1;
			
			//Inc byte counters
			circBufPtr = (circBufPtr + 1) % m;
			nChar += 1;
			if (nChar >= 4) {
				bInd += 1;
				nChar = nChar % 4;
			}
		}
		return map;
	}
	
	
	
	//Get The next non-zero posn in mapStartPosn and then get its value
	private static int getNextStart(int posn) {
		int res = hashCount;
		for (int i =posn+1;i< mapStartPosn.length;i++)
			
			if (Integer.compareUnsigned(mapStartPosn[i],0) > 0) {
				res = mapStartPosn[i];
				break;
			}
		return res;
	}
	
	//num = kmer
	public static int stringPopularity(byte[] chrA,int seqStPos, int m, int[] splitA,int num) {
		byte[] pat = br.getBytesFromBB(chrA,seqStPos,m);
		int totalHits = 0;

		int mapStart = 0;
		int mapNext = 0;	
		int hp = 0;
		
		for (int startPosn: splitA) {
			hp = expHash(pat,startPosn,num);

			//curr string hashes to an empty bucket ie. no hits found for this string
			if ((hp != firstMapPosn) && (mapStartPosn[hp] == 0)) {
				mapStart = mapNext = 0;
			}else {	
				mapStart = mapStartPosn[hp];
				mapNext = getNextStart(hp);
			}
			totalHits += (mapNext - mapStart);
		}
		return totalHits;
	}

	public static TreeMap<Integer,int[]> querySplit(byte[] pat,int[] splitA,int num) {
			
		int mapStart = 0;
		int mapNext = 0;
		int hp = 0;
		
		//Get Probables for each subset
		TreeMap<Integer,int[]> probables = new TreeMap<Integer,int[]> ();
		for (int startPosn: splitA) {
			hp = expHash(pat,startPosn,num);
			//curr string hashes to an empty bucket ie. no hits found for this string
			if ((hp != firstMapPosn) && (mapStartPosn[hp] == 0)) {
				mapStart = mapNext = 0;
			}else {	
				mapStart = mapStartPosn[hp];
				mapNext = getNextStart(hp);
			}
			int hits = mapNext - mapStart;
			
			//To restrict number of hits and spped up alignment
			//if (hits > 256)
				//continue;
			if (dbg > 0)
				System.out.println("HashIndex:querySplit() hp=" + hp + " mapStart=" + mapStart +" mapNext=" + mapNext  + " startPosn=" + startPosn + " hits=" + hits);
		
			long tt = Integer.toUnsignedLong(mapStart * Integer.BYTES);
			tt += mapValues;
			
			for (int i=0;Integer.compareUnsigned(i,hits)<0;i++) {
				
				int hit = unsafe.getInt(tt) - startPosn;
				tt += Integer.BYTES;
				//For DBG detail
				//if (dbg > 1)
					//System.out.println("DBG Detail mapV=" + mapValues[mapStart+i] + " adjusted=" + (mapValues[mapStart+i] - startPosn) + " startPosn=" + startPosn);
				//Add to probables
				int[] val= probables.get(hit);
				if (val == null) 
					probables.put(hit,new int[]{1});
				else
					++val[0];
			}
		}

		//SKIP check if these are true hits
		/**
        for pr in prob:
            if p == tA[pr:pr+m]:
                hits.append(pr)
		**/
        return probables;
	}
	
	public static void printFreqProb(int[] fProb) {
		System.out.println("Total Entries=" + IntStream.of(fProb).parallel().sum());
		for (int x : fProb)
			System.out.print(x + " ");
		System.out.println();
	}
      
	public static void printMap(int[] map) {
		for (int i=0; i< map.length; i++) {
			System.out.println("i=" + i + " count=" + map[i]);
		}
	}
	
	public static void printCounter(Map<Integer,int[]> ctr) {
		System.out.println("printCounter len keys = " + ctr.keySet().size());
		for (int x : ctr.keySet()) {
			System.out.println("key= " + x + " val=" + ctr.get(x)[0]);
		}
	
	}
	
}