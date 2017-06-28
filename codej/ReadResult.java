/**
 * ReadResult: result of a read alignment
 */
public class ReadResult {
	public double confidence;
	public int r;
	
	ReadResult(double conf,int r) {
		this.confidence = conf;
		this.r = r;
	}
};
