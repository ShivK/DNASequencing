import java.nio.charset.StandardCharsets;
import java.util.Arrays;


public class BinaryRep {
	private static final int[] g_codes = {65,67,71,84};
	private static final byte[] g_vals = {0b00000000,0b00000001,0b00000010,0b00000011};
	
	private static final byte[] complement = {0b00000011,0b00000010,0b00000001,0b00000000};
	
	public static final byte[] g_masks = {(byte)0b11000000,0b00110000,0b00001100,0b00000011};

	//get kmer charcters from the ith char from Binary packed Array (4 chars to a byte)
	//the kmer chars will be treturned in a byte array
	//Assume that i < lenbA ie the number of chars in bA
	public static byte[] getBytesFromBB(byte[] bA,int i,int kmer) {
		byte[] retA = new byte[kmer];
		int retInd = 0;

		//get the start byte point
		int bInd = Integer.divideUnsigned(i,4);
		//get the ith char in the mask
		int n = Integer.remainderUnsigned(i,4);
		//get the byte
		byte b = bA[bInd];
    
		while (retInd < kmer) {
			//get the byte
			//#right shift
			byte rb = (byte)(((b & g_masks[n])&0xff) >> ((3 - n) * 2));
			/**
			System.out.println("retInd=" + retInd);		
			System.out.println(Integer.toBinaryString((b & 0xFF) + 0x100).substring(1));
			System.out.println("n=" +n);		
			System.out.println(Integer.toBinaryString((rb & 0xFF) + 0x100).substring(1));
			**/
			retA[retInd] = rb;
			retInd += 1;
			n += 1;
			if ((n % 4) == 0) {
				bInd += 1;
				n = 0;
				b = bA[bInd];
			}
		}
		return retA;
	}

	/**
	*convert a line string of genomic codes to binary code - 1 per byte
	* input cd - line of genome string
	* Assume All valid genomic codes
	* Assume length of bA is equal to len(cd)
	**/
	public byte[] convGLineToBBytes(String cd,byte[] bA) {
		
		for (int i =0; i< cd.length(); i++) {
			int c = (int)cd.charAt(i);
			if (c == 65)  //A
				bA[i] = g_vals[0];
			 else if (c == 67) 
				bA[i] = g_vals[1];
			else if (c == 71) 
				bA[i] = g_vals[2];
			else if (c == 84)
				bA[i] = g_vals[3];
			else 
				bA[i] = g_vals[0];
		}
		return bA;
	}

	/**
	#s - sequence
	# t - byteaeeat created for the reversecomplement
	# i = len(t) -1
	**/
	byte[] reverseComplementBBytes(byte[]src,byte[] dst) {
		int num = dst.length - 1;
		for (byte base : src) {
			dst[num] = complement[base];
			num-=1;
		}
		return dst;
	}

	/**
	#convert a line string of genomic codes to binary
	# input cd - line of genome string
	# if genome is N, map it to A and note the posn
	* currently noy noting the N posn
	**/
	public byte[] convGLineToB(String cd,byte[] bA,int stB,int[] NCtA) {
		//System.out.println("convGlineToB:" + cd + " stB=" + stB);
		//printBin(bA);

		//NStCtr = 0
		//NEndCtr = 0
		//NCount = False
		int btCtr = 0;
		int limit = 0;
		for (int i =0; i< cd.length(); i+= 4) {
			//wd = cd[i:i+4]
			int shift = 6;
			limit =  cd.length() - i;
			if (limit > 4)
				limit = 4;
			
			//System.out.println("convGlineToB: limit=" + limit);
			byte setB;
			for (int j=0; j< limit ; j++) {
				setB = 0b00000000;
				
				int c = (int)cd.charAt(i+j);
				
				if (c == 65) { //A
					setB = (byte) (g_vals[0]  << shift);
				} else if (c == 67) 
					setB = (byte) (g_vals[1]  << shift);
				else if (c == 71) 
					setB = (byte) (g_vals[2]  << shift);
				else if (c == 84)
					setB = (byte) (g_vals[3]  << shift);
				else {
					setB = (byte) (g_vals[0]  << shift);
					//set N
				}
				bA[stB+btCtr] |= setB;
				shift -= 2;
				
			}
					
			btCtr += 1;
		}


		return bA; //,NCtA
	}

	public static void printBin(byte[] bA) {
		for (byte b : bA) {
			System.out.println(Integer.toBinaryString((b & 0xFF) + 0x100).substring(1));
		}
	}
	
	//bA is in B format, st - start number in bA, kmer - k chars from i
	public static void printBString(byte[] bA, int st,int kmer) {
		byte[] charA = getBytesFromBB(bA,st,kmer);
		printString(charA);
	}

	public static String getString(byte[] charA,int ln) {
		StringBuffer sb = new StringBuffer(ln);
		for (int i=0;i<ln;i++) 
			sb.append((char)g_codes[charA[i]]);
		return sb.toString();
	}
	
	public static void printString(byte[] charA) {
		for (int i=0;i<charA.length;i++) 
			System.out.print((char)g_codes[charA[i]]);
		System.out.println();
	}
		
	
	public static byte[] hexStringToByteArray(String s) {
		int len = s.length();
		byte[] data = new byte[len / 2];
		System.out.println("hex=" + s);
		
		for (int i = 0; i < len; i += 1) {
			System.out.println(s.charAt(i));
			System.out.println((int)s.charAt(i));
		}
		return data;
	}
	
	public byte[] getBytes(String seq) {
		byte[] b = seq.getBytes(StandardCharsets.US_ASCII); 
		System.out.println("chrSeq=" + seq);
		System.out.println("chrSeq bytes=" + Arrays.toString(b) + "bL=" + b.length);
		String string = new String(b, StandardCharsets.US_ASCII);
		System.out.println("chrSeq bytst=" + string);
		//Binary String
		printBin(b);
		
		
		return b;
	}
	

      
      
}