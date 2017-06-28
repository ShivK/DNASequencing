/**
 * Position: describe the probable position of a read within the genome
 */
import java.text.DecimalFormat;

public class Probable {
	public int chrId = 20;
	public int diff;
	public int reada; //position of read strand 1
	public int readb; //position of read strand 2
	public int countA;
	public int countB;
	public int total;    //countA + count B + tolerance count
	public int tolCount; // entries within tolerance count
	public int type; // 0 -> +/- 1-> -/+
	public float conf;

	Probable() {
	}
	
	void printDebug() {
		System.out.println("chrId = " + chrId + " reada=" + reada + " diff=" + diff + " countA=" + countA + " countB=" + countB + " tolCount=" + tolCount + " total=" + total + " conf=" + conf);
	}
	
	void setValues(int diff,int reada,int readb,int countA,int countB, int total,int type) {
		this.diff= diff;
		this.reada = reada;
		this.readb = readb;
		this.countA = countA;
		this.countB = countB;
		this.total = total;
		this.type = type;
	}

	void setRealChr(int chrId,int reada,int readb) {
		this.chrId= chrId;
		this.reada = reada;
		this.readb = readb;
	}
	
	//Type == 0 -> M1 1 -> M2
	String[] getResStr() {
		DecimalFormat df = new DecimalFormat("#.##");

		String[] res = new String[2];
		
		StringBuffer sb1 = new StringBuffer();
		StringBuffer sb2 = new StringBuffer();
		sb1.append(",");
		sb1.append(chrId);
		sb1.append(",");
		sb1.append(reada+1);
		sb1.append(",");
		sb1.append(reada + 150);
		sb1.append(",");

		sb2.append(",");
		sb2.append(chrId);
		sb2.append(",");
		sb2.append(readb+1);
		sb2.append(",");
		sb2.append(readb + 150);
		sb2.append(",");
		
		if (type == 0) {
			sb1.append("+");
			sb2.append("-");
		}else {
			sb1.append("-");
			sb2.append("+");
		}
		sb1.append(",");
		sb2.append(",");
		
		sb1.append(df.format(conf));
		sb2.append(df.format(conf));
		
		res[0] = sb1.toString();
		res[1] = sb2.toString();
		
		return res;

	}
	
};
