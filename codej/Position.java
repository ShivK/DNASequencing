/**
 * Position: describe the position of a read within the genome
 */


public class Position {
	public int chrName; //eg: chr - 20
	public int from;
	public int to;
	public char strand;
	
	Position(int rn,int from, int to,char st) {
		this.chrName = rn;
		this.from = from;
		this.to = to;
		this.strand = st;
	}
	
};
