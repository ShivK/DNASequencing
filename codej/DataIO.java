import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

import java.nio.file.Files;
import java.io.IOException;
import java.nio.file.Paths;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;

import java.io.RandomAccessFile;
import java.nio.channels.FileChannel;





public class DataIO {
	
	//Read the full file
	private String[] readFileInt(String fname, boolean fasta) {
		try {
			ArrayList<String> fileA = new ArrayList<String>();
			FileInputStream fin = null;
			BufferedReader reader = null;
			
			fin = new FileInputStream(fname);
			reader = new BufferedReader(new InputStreamReader(fin));

			String line = null;
			if (fasta)
				line = reader.readLine(); //skip header

			
			line = reader.readLine(); 
			while(line != null){
				fileA.add(line);
				line = reader.readLine();
			}           
			reader.close();
			fin.close();
			return fileA.toArray(new String[0]);

		}catch(Exception e) {
			System.out.println(e);
			return null;
		}

	}

	public String[] readFile(String fname) {
		return readFileInt(fname,false);
	}

	public String[] readFastaFile(String fname) {
		return readFileInt(fname,true);
	}
	
	public int readSeqFile(String fa1Path,String fa2Path,ArrayList<String> readName,ArrayList<String> readSequence) {
		try {
			//reading file line by line using BufferedReader       
			FileInputStream fin1 = new FileInputStream(fa1Path);
			BufferedReader reader1 = new BufferedReader(new InputStreamReader(fin1));

			FileInputStream fin2 = new FileInputStream(fa2Path);
			BufferedReader reader2 = new BufferedReader(new InputStreamReader(fin2));
			
			String line1 = reader1.readLine(); //name line
			String line2 = reader2.readLine(); //name line
			
			while(line1 != null && line2 != null){
				//Add name
				readName.add(line1.substring(1,line1.length()));
				readName.add(line2.substring(1,line2.length()));
				//System.out.println(line);
				line1 = reader1.readLine(); //read seq
				line2 = reader2.readLine();
				readSequence.add(line1);
				readSequence.add(line2);
				
				line1 = reader1.readLine(); //read name 
				line2 = reader2.readLine();
			}           
			reader1.close();
			reader2.close();
			fin1.close();
			fin2.close();
			return 0;

		}catch(Exception e) {
			System.out.println(e);
			return -1;
		}

	}
	
	//static ByteBuffer buffer = ByteBuffer.allocate(813920000); 
 	public static byte[] loadBinaryFile(String filename,long start,int length) {
		try { 
			// Opening RandomAccessFile for reading data 
			RandomAccessFile store = new RandomAccessFile(filename, "r"); 
			// getting file channel 
			FileChannel channel = store.getChannel(); 
			// preparing buffer to read data from file 
			ByteBuffer buffer = ByteBuffer.allocate(length); 
			// reading data from file channel into buffer 
			int numOfBytesRead = channel.read(buffer,start); 
			//System.out.println("number of bytes read : " + numOfBytesRead); 
			// You need to filp the byte buffer before reading 
			buffer.flip(); 
			
			channel.close(); 
			store.close(); 
			return buffer.array();
		} catch (IOException e) { 
			e.printStackTrace(); 
			return new byte[0];
		} 
		
	}
	
	
	public static int saveBinaryFile(String filename, byte[] data) {
	  try {
		Files.write(Paths.get(filename), data);
		return 0;
	  } catch (IOException e) {
		  System.out.println("Unable to write file: " + filename + "\n");
		return -1;
	  }
	}
	
	public static int saveIntFile(String filename, int[] data) {  
		return saveBinaryFile(filename, toByteArrayInt(data));
	}	

	private static byte[] toByteArrayInt(int[] data) {
		ByteBuffer byteBuffer = ByteBuffer.allocate(data.length * 4);        
        IntBuffer intBuffer = byteBuffer.asIntBuffer();
        intBuffer.put(data);
		System.out.println("toByteArrayInt: array size = " + byteBuffer.array().length + " capacity = " + byteBuffer.capacity());

        return byteBuffer.array();	
	}
	
	
	private static int[] toIntArray(byte[] byteArray){
	  int times = Integer.SIZE / Byte.SIZE;
	  
	  int[] ints = new int[byteArray.length / times];
	  for (int i = 0; i < ints.length; i++){
		ints[i] = ByteBuffer.wrap(byteArray, i * times, times).getInt();
	  }
	  return ints;
	}
	
	
	public static int[] loadFromIntFile(String filename,long start,int length) {
		if (start ==0 && length == -1)
			return toIntArray(loadBinaryFile(filename,start,1024));
		else {
			return toIntArray(loadBinaryFile(filename,start*4,length*4));
		}
	}


    public static void main(String args[]) {
		System.out.println("dataIO.start()" + args.length);
		String fname = args[0];
		DataIO d = new DataIO();
		d.readFile(fname);
	}
      
      
}